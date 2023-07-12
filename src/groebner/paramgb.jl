
"""
    paramgb(blackbox) -> basis

## Example

```jldoctest
using ParamPunPam, Nemo

Rparam, (a, b) = PolynomialRing(QQ, ["a", "b"])
R, (x, y, z) = PolynomialRing(FractionField(Rparam), ["x", "y", "z"], ordering=:degrevlex)

F = [
    x^2 + x + (a + 1),
    x*y + b*y*z + 1//(a*b),
    x*z + z + b
]

paramgb(BasicBlackboxIdeal(F))
```

"""
function paramgb(
        blackbox::T;
        kwargs...
        ) where {T<:AbstractBlackboxIdeal}
    up_to_degree = get(kwargs, :up_to_degree, (div(typemax(Int), 2), div(typemax(Int), 2)))
    guess_degrees = get(kwargs, :guess_degrees, :no)
    rational_interpolator = get(kwargs, :rational_interpolator, VanDerHoevenLecerf())
    polynomial_interpolator = get(kwargs, :polynomial_interpolator, :primes_benortiwari)
    @info "Computing parametric GB" up_to_degree rational_interpolator polynomial_interpolator    
    _paramgb(
        blackbox,
        up_to_degree,
        guess_degrees, rational_interpolator, polynomial_interpolator
    )
end

function _paramgb(
        blackbox,
        up_to_degree,
        guess_degrees, rational_interpolator, polynomial_interpolator
    )
    # The struct to keep track of modular computation related stuff
    modular = ModularTracker(blackbox)
    # The struct to store the state of the computation
    state = GroebnerState(blackbox)
    # Discover the shape of the groebner basis:
    # its size, and the sizes of polynomials it contains
    discover_shape!(state, modular, η=2)
    ### guess_degrees
    discover_param_degrees!(state, modular)
    # Interpolate the exponents in the parametric coefficients
    # (this uses exactly 1 prime number)
    ### rational_interpolator, polynomial_interpolator
    interpolate_param_exponents!(
        state, modular, up_to_degree, rational_interpolator)
    # Interpolate the rational coefficients of the parametric coefficients
    # (this is expected to use 1 prime number, but may use more)
    recover_coefficients!(state, modular)
    # Combine and return the above two 
    basis = construct_basis(state)
    basis
end

@noinline __throw_unlucky_cancellation() = throw(AssertionError("Unlucky cancellation of coefficients!"))

# Discovers shape of the groebner basis of the ideal from `state`
# by specializing it at a random point (preferably, modulo a prime).
#
# If η > 0 is given, the algorithm will confirm the shape
# by additionaly specializing it at η random points.
function discover_shape!(state, modular; η=2)
    iszero(η) && (@warn "Fixing the shape of the basis from 1 point is adventurous.")
    @info "Specializing at $(1) + $(η) random points to guess the basis shape.."
    blackbox = state.blackbox
    # Guess the shape for 1 lucky prime:
    reduce_mod_p!(blackbox, modular.ff)
    @label Start
    # specialize at a random lucky point and compute GBs
    randompoints = map(_ -> randluckyspecpoint(state, modular.ff), 1:1 + η)
    polysspecmodp = map(point -> evaluate_mod_p(blackbox, point), randompoints)
    graph, gb = groebner_learn(polysspecmodp[1])
    state.graph = graph
    bases = empty(polysspecmodp)
    for F in polysspecmodp
        flag, gb = groebner_apply!(graph, F)
        if !flag
            __throw_unlucky_cancellation()
        end
        push!(bases, gb)
    end
    # decide the "right" basis according to the major rule
    success, basis = majorrule(bases)
    # if not all of the bases have the same shape
    if !success
        η = 2*(1 + η)
        @goto Start
    end 
    state.shape = basisshape(first(bases))
    @info "The shape of the basis is: $(length(basis)) polynomials with monomials" state.shape
    nothing
end

is_interpolated_heuristic(dD, DD, dN, DN) = (dD + 3) < div(3DD, 4) && (dN + 3) < div(3DN, 4)

function discover_param_degrees!(state, modular)
    @info "Specializing at random points to guess the total degrees in parameters.."
    blackbox = state.blackbox
    Ru, _ = PolynomialRing(modular.ff, :u)
    K = base_ring(Ru)
    Rx = parent(blackbox)
    Ra = base_ring(Rx)
    n = length(gens(Ra))
    shift = distinct_points(K, n)
    shape = state.shape
    graph = state.graph

    N, D = 1, 1
    npoints = N + D + 2

    coeffs = Vector{Vector{Vector{elem_type(K)}}}(undef, length(shape))
    interpolated = Vector{Vector{Tuple{elem_type(Ru), elem_type(Ru)}}}(undef, length(shape))
    degrees = Vector{Vector{Tuple{Int, Int}}}(undef, length(shape))
    for i in 1:length(shape)
        interpolated[i] = Vector{Tuple{elem_type(Ru), elem_type(Ru)}}(undef, length(shape[i]))
        degrees[i] = Vector{Tuple{Int, Int}}(undef, length(shape[i]))
        coeffs[i] = Vector{Vector{elem_type(K)}}(undef, length(shape[i]))
        for j in 1:length(shape[i])
            coeffs[i][j] = Vector{elem_type(K)}(undef, npoints)
        end
    end

    J = 1
    univ_x_points = Vector{elem_type(K)}()
    all_interpolated = false
    while !all_interpolated
        all_interpolated = true
        N, D = 2N, 2D
        interpolator = FasterCauchy(Ru, N, D)
        npoints = N + D + 2
        @info "Using $npoints points.."
        univ_x_points = distinct_points(K, npoints - length(univ_x_points), univ_x_points)
        x_points = map(point -> repeat([point], n) .+ shift, univ_x_points)
        for i in 1:length(shape)
            for j in 1:length(shape[i])
                resize!(coeffs[i][j], npoints)
            end
        end
        for idx in J:npoints
            point = x_points[idx]
            Ip = evaluate_mod_p(blackbox, point)
            flag, basis = groebner_apply!(graph, Ip)
            if !flag
                __throw_unlucky_cancellation()
            end
            # this to some extent duplicates the code from groebner_apply!
            !check_shape(shape, basis) && @warn "Bad substitution!"
            
            for i in 1:length(coeffs)
                for j in 1:length(coeffs[i])
                    coeffs[i][j][idx] = coeff(basis[i], j)
                end
            end
        end
        for i in 1:length(coeffs)
            for j in 1:length(coeffs[i])
                P, Q = interpolate!(interpolator, univ_x_points, coeffs[i][j])
                interpolated[i][j] = (P, Q)
                degrees[i][j] = (degree(P), degree(Q))
                dp, dq = degree(P), degree(Q)
                if !is_interpolated_heuristic(dp, N, dq, D)
                    all_interpolated = false
                end
            end
        end
        J = npoints + 1
    end
    state.param_degrees = degrees
    @info "Success! $(npoints) points used."
    @info "The total degrees in the coefficients" state.param_degrees
    nothing
end

function get_interpolator(interpolator::CuytLee, Ru, Nd, Dd, Nds, Dds, Nt, Dt)
    FasterCuytLee(
        Ru, Nd, Dd, Nds, Dds, Nt, Dt)
end

function get_interpolator(interpolator::VanDerHoevenLecerf, Ru, Nd, Dd, Nds, Dds, Nt, Dt)
    FasterVanDerHoevenLecerf(
        Ru, Nd, Dd, Nds, Dds, Nt, Dt)
end

function interpolate_param_exponents!(
        state, modular, up_to_degree, rational_interpolator)
    @info "Interpolating the exponents in parameters.."
    blackbox = state.blackbox
    Rx = parent(blackbox)
    Ra = base_ring(Rx)
    Ru, _ = PolynomialRing(modular.ff, symbols(Ra), ordering=Nemo.ordering(Ra))
    K = base_ring(Ru)
    n = length(gens(Ra))
    shape = state.shape
    graph = state.graph

    Nt, Dt = 1, 1
    total_degrees = state.param_degrees
    Nd0 = maximum(d -> maximum(dd -> dd[1], d), total_degrees)
    Dd0 = maximum(d -> maximum(dd -> dd[2], d), total_degrees)
    
    Nd = min(Nd0, up_to_degree[1] + 1)
    Dd = min(Dd0, up_to_degree[2] + 1)
    
    Nds, Dds = repeat([Nd], n), repeat([Dd], n)
    npoints = 0

    param_exponents = Vector{Vector{Tuple{elem_type(Ru), elem_type(Ru)}}}(undef, length(shape))
    coeffs = Vector{Vector{Vector{elem_type(K)}}}(undef, length(shape))
    for i in 1:length(shape)
        param_exponents[i] = Vector{Tuple{elem_type(Ru), elem_type(Ru)}}(undef, length(shape[i]))
        coeffs[i] = Vector{Vector{elem_type(K)}}(undef, length(shape[i]))
        for j in 1:length(state.shape[i])
            coeffs[i][j] = Vector{elem_type(K)}(undef, npoints)
        end
    end

    interpolator = get_interpolator(rational_interpolator, Ru, Nd, Dd, Nds, Dds, Nt, Dt)

    J = 1
    all_interpolated = false
    @info "Interpolating for degrees:\nnumerator $Nd, denominator $Dd"
    while !all_interpolated
        all_interpolated = true

        x_points = get_evaluation_points!(interpolator)
        npoints = length(x_points)
        @info "Using $npoints points.."
        
        for i in 1:length(shape)
            for j in 1:length(shape[i])
                resize!(coeffs[i][j], npoints)
            end
        end

        for idx in J:npoints
            point = x_points[idx]
            Ip = evaluate_mod_p(blackbox, point)
            
            flag, basis = groebner_apply!(graph, Ip)
            if !flag
                __throw_unlucky_cancellation()
            end
            !check_shape(shape, basis) && @warn "Bad substitution!"

            # TODO: to be removed
            # @assert basisshape(basis) == shape

            for i in 1:length(shape)
                for j in 1:length(shape[i])
                    coeffs[i][j][idx] = coeff(basis[i], j)
                end
            end
        end
        
        for i in 1:length(shape)
            for j in 1:length(shape[i])
                if total_degrees[i][j][1] > up_to_degree[1] || total_degrees[i][j][2] > up_to_degree[2]
                    param_exponents[i][j] = (Ru(-1), Ru(-1))
                    continue
                end
                P, Q = interpolate!(interpolator, coeffs[i][j])
                param_exponents[i][j] = (P, Q)
                dp, dq = total_degree(P), total_degree(Q)
                if !(dp >= total_degrees[i][j][1] && dq >= total_degrees[i][j][2])
                    all_interpolated = false
                    break
                end
            end
        end
        J = npoints + 1
    end
    state.param_exponents = param_exponents
    @info "Success! $(npoints) points used."
    maxDn = maximum(l -> maximum(total_degree ∘ first, l), param_exponents)
    maxDd = maximum(l -> maximum(total_degree ∘ last, l), param_exponents)
    maxTn = maximum(l -> maximum(length ∘ first, l), param_exponents)
    maxTd = maximum(l -> maximum(length ∘ last, l), param_exponents)
    @info "Output summary:
    Maximal interpolated degrees are: $maxDn for num. and $maxDd for den.
    Maximal number of interpolated terms are: $maxTn for num. and $maxTd for den.
    "
    nothing
end

function recover_coefficients!(state, modular)
    @info "Recovering the coefficients.."
    blackbox = state.blackbox
    Rorig = parent(blackbox)
    Rparam = parent_params(blackbox)
    Rorig_frac, _ = PolynomialRing(Nemo.FractionField(Rparam), symbols(Rorig), ordering=Nemo.ordering(Rorig))
    Rparam_frac = base_ring(Rorig_frac)
    polysreconstructed = Vector{elem_type(Rorig_frac)}(undef, length(state.shape))
    p = convert(Int, characteristic(modular.ff))
    for i in 1:length(state.shape)
        coeffsrec = Vector{elem_type(Rparam_frac)}(undef, length(state.shape[i]))
        for j in 1:length(state.shape[i])
            P, Q = state.param_exponents[i][j]
            Prec = map_coefficients(c -> rational_reconstruction(Int(data(c)), p), P)
            Qrec = map_coefficients(c -> rational_reconstruction(Int(data(c)), p), Q)
            coeffsrec[j] = Prec // Qrec
        end
        polysreconstructed[i] = Rorig_frac(coeffsrec, map(e -> exponent_vector(e, 1), state.shape[i]))
    end
    state.param_coeffs = polysreconstructed
    @info "Success! Used $(1) prime in total :)"
    nothing
end

function construct_basis(state)
    state.param_coeffs
end

function check_shape(best, current)
    if length(best) != length(current)
        return false
    end
    if map(length, best) != map(length, current)
        return false
    end
    # TODO: remove this check
    # if best != map(f -> collect(monomials(f)), current)
    #     return false
    # end
    true
end

function majorrule(bases::Vector{T}) where {T}
    if length(bases) == 1
        return true, first(bases)
    end
    lengthes = map(length, bases)
    if length(unique!(lengthes)) > 1
        @warn "One of the computed bases has a different shape"
        return false, first(bases)
    end
    true, first(bases)
end

function basisshape(basis)
    map(collect ∘ monomials, basis)
end
