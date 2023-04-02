
"""
    paramgb

Given an array of polynomials `polys` over a field of rational functions
computes the Groebner basis of the ideal generated by `polys`.

The algorithm is probabilistic and succeeds with a high probability.

Examples:

```julia
using Nemo, ParamPunPam

Rparam, (a, b) = PolynomialRing(QQ, ["a", "b"], ordering=:degrevlex)
R, (x, y, z) = PolynomialRing(FractionField(Rparam), ["x", "y", "z"], ordering=:degrevlex)

ParamPunPam.paramgb([a*x^2 + 1, y^2*z + (1//b)*y])
```

Supported keyword arguments:

- `guess_degrees`, try to guess the degrees of parameters in the output.
Possible options are either `:total_degrees`, or `:partial_degrees`.
- `rational_interpolator`, the rational function interpolation algorithm to use.
Possible options are `:cuytlee` and `:vdhoevenlecerf`.
- `polynomial_interpolator`, the polynomial interpolation algorithm to use.
Possible options are `:primes_benortiwari` and `:kron_benortiwari`.

"""
function paramgb(
        polys::Vector{T};
        kwargs...
        ) where {T}
    guess_degrees = get(kwargs, :guess_degrees, :no)
    rational_interpolator = get(kwargs, :rational_interpolator, :cuytlee)
    polynomial_interpolator = get(kwargs, :polynomial_interpolator, :primes_benortiwari)    
    metainfo = peek_at_input(polys)
    unified_polys = unify_input(polys, metainfo)
    _paramgb(
        unified_polys, metainfo,
        guess_degrees, rational_interpolator, polynomial_interpolator
    )
end

# Returns the ring of polynomials in parameters
getparamring(coeffring) = throw(DomainError(coeffring, "Unknown coefficient ring."))
getparamring(coeffring::Nemo.FracField) = true, base_ring(coeffring)
getparamring(coeffring::Nemo.MPolyRing) = false, coeffring

# Checks the input and returns some meta-information about it
function peek_at_input(polys)
    @assert !isempty(polys) "Empty input is invalid."
    Rx = parent(first(polys))
    over_fractions, Rparam = getparamring(base_ring(Rx))
    K = base_ring(Rparam)
    @assert all(x -> parent(x) == Rx, polys) "All polynomials must be in the same ring."
    @assert typeof(K) === Nemo.FlintRationalField || typeof(K) === typeof(AbstractAlgebra.QQ) "Coefficient ring must be Nemo.QQ or AbstractAlgebra.QQ"
    (over_fractions=over_fractions,)
end

#=
    Does something cool with the input.
=#
function unify_input(polys, metainfo)
    filter!(!iszero, polys)
end

function _paramgb(
        polys, metainfo,
        guess_degrees, rational_interpolator, polynomial_interpolator
    )
    # The struct to keep track of modular computation related stuff
    modular = ModularTracker(polys)
    # The struct to store the state of the computation
    state = GroebnerState(polys)
    # Discover the shape of the groebner basis:
    # its size, and the sizes of polynomials it contains
    discover_shape!(state, modular, η=2)
    ### guess_degrees
    discover_param_degrees!(state, modular)
    # Interpolate the exponents in the parametric coefficients
    # (this uses exactly 1 prime number)
    ### rational_interpolator, polynomial_interpolator
    interpolate_param_exponents!(state, modular, rational_interpolator)
    # Interpolate the rational coefficients of the parametric coefficients
    # (this is expected to use 1 prime number, but may use more)
    recover_coefficients!(state, modular)
    # Combine and return the above two 
    basis = construct_basis(state)
    basis
end

# Discovers shape of the groebner basis of the ideal from `state`
# by specializing it at a random point (preferably, modulo a prime).
#
# If η > 0 is given, the algorithm will confirm the shape
# by additionaly specializing it at η random points.
function discover_shape!(state, modular; η=2)
    iszero(η) && (@warn "Fixing the shape of the basis from 1 point is adventurous.")
    @info "Specializing at $(1) + $(η) random points to guess the basis shape.."
    # Guess the shape for 1 lucky prime:
    polysmodp = reducemodp(state.polys_fracfree, modular)
    @label Start
    # specialize at a random lucky point and compute GBs
    randompoints = map(_ -> randluckyspecpoint(state, modular.ff), 1:1 + η)
    polysspecmodp = map(point -> specialize(polysmodp, point), randompoints)
    @assert all(F -> ordering(parent(first(F))) === :degrevlex, polysspecmodp)
    bases = map(F -> groebner(F, linalg=:prob), polysspecmodp)
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

is_interpolated_heuristic(dD, DD, dN, DN) = (dD + 3) < div(2DD, 4) && (dN + 3) < div(2DN, 4)

function discover_param_degrees!(state, modular)
    @info "Specializing at random points to guess the total degrees in parameters.."
    Ru, _ = PolynomialRing(modular.ff, :u)
    K = base_ring(Ru)
    Rx = parent(first(state.polys_fracfree))
    Ra = base_ring(Rx)
    n = length(gens(Ra))
    polysmodp = reducemodp(state.polys_fracfree, modular)
    shift = distinct_points(K, n)
    shape = state.shape

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
            Ip = specialize(polysmodp, point)
            @assert ordering(parent(first(Ip))) === :degrevlex
            basis = groebner(Ip, linalg=:prob)
            
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

function interpolate_param_exponents!(
        state, modular, rational_interpolator)
    @info "Interpolating the exponents in parameters.."
    Rx = parent(first(state.polys_fracfree))
    Ra = base_ring(Rx)
    Ru, _ = PolynomialRing(modular.ff, symbols(Ra))
    K = base_ring(Ru)
    n = length(gens(Ra))
    shape = state.shape
    polysmodp = reducemodp(state.polys_fracfree, modular)

    Nt, Dt = 1, 1
    total_degrees = state.param_degrees
    Nd = maximum(d -> maximum(dd -> dd[1], d), total_degrees)
    Dd = maximum(d -> maximum(dd -> dd[2], d), total_degrees)
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

    interpolator = FasterVanDerHoevenLecerf(
        Ru, Nd, Dd, Nds, Dds, Nt, Dt)

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
            Ip = specialize(polysmodp, point)
            @assert ordering(parent(first(Ip))) === :degrevlex
            
            basis = groebner(Ip, linalg=:prob)
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
                P, Q = interpolate!(interpolator, coeffs[i][j])
                param_exponents[i][j] = (P, Q)
                dp, dq = total_degree(P), total_degree(Q)
                if !(dp >= total_degrees[i][j][1] && dq >= total_degrees[i][j][2])
                    all_interpolated = false
                end
            end
        end
        J = npoints + 1
    end
    state.param_exponents = param_exponents
    @info "Success! $(npoints) points used."
    @info "The exponents in the coefficients" state.param_exponents
    nothing
end

function recover_coefficients!(state, modular)
    @info "Recovering the coefficients.."
    Rorig = parent(first(state.polys))
    Rparam = base_ring(Rorig)
    Ra = base_ring(Rx)
    n = length(gens(Ra))
    polysreconstructed = Vector{elem_type(Rorig)}(undef, length(state.shape))
    p = convert(Int, characteristic(modular.ff))
    for i in 1:length(state.shape)
        coeffsrec = Vector{elem_type(Rparam)}(undef, length(state.shape[i]))
        for j in 1:length(state.shape[i])
            P, Q = state.param_exponents[i][j]
            Prec = map_coefficients(c -> rational_reconstruction(Int(data(c)), p), P)
            Qrec = map_coefficients(c -> rational_reconstruction(Int(data(c)), p), Q)
            coeffsrec[j] = Prec // Qrec
        end
        polysreconstructed[i] = Rorig(coeffsrec, map(e -> exponent_vector(e, 1), state.shape[i]))
    end
    state.param_coeffs = polysreconstructed
    @info "Success! Used $(1) prime in total :)"
    nothing
end

function reducemodp(polys, modular::ModularTracker)
    ff = modular.ff
    @info "Reducing modulo $(ff).."
    polysmodp = map(
        poly -> map_coefficients(
            f -> map_coefficients(
                c -> ff(c), 
                f
            ), 
            poly
        ), 
        polys
    )
    polysmodp
end

function construct_basis(state)
    state.param_coeffs
end

function specialize(polys, point)
    map(f -> map_coefficients(c -> evaluate(c, point), f), polys)
end

function check_shape(best, current)
    if length(best) != length(current)
        return false
    end
    if map(length, best) != map(length, current)
        return false
    end
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
