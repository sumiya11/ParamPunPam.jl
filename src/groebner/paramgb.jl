
"""
    paramgb(polys; options...) -> basis

Computes the Groebner basis of the ideal generated by `polys`, where `polys` are
polynomials over a field of rational functions.

The algorithm is probabilistic and succeeds with a high probability.

## Possible Options

Supported keyword arguments:

- `ordering`: monomial ordering. Same as for Groebner.jl, for example, `Lex()` or
  `DegRevLex()`.
- `up_to_degree`: Compute parameters up to a fixed total degree of
  numerators and denominators e.g., `up_to_degree=(2, 2)`.
- `rational_interpolator`: The rational function interpolation algorithm.
  Possible options are `:CuytLee` and `:VanDerHoevenLecerf` (default).
- `estimate_degrees`: If `true`, estimates the total degrees of parameters
  before starting the interpolation. Default is `true`.
- `assess_correctness`: If `true`, check that the basis is correct with high
  probability. Default is `false`.

## Example

```jldoctest
using Nemo, ParamPunPam

Rparam, (a, b) = PolynomialRing(QQ, ["a", "b"])
R, (x, y, z) = PolynomialRing(FractionField(Rparam), ["x", "y", "z"], ordering=:degrevlex)

ParamPunPam.paramgb([a*x^2 + 1, y^2*z + (1//b)*y])
```
"""
function paramgb(polys::Vector{T}; kwargs...) where {T}
    paramgb(BasicBlackboxIdeal(polys); kwargs...)
end

"""
    paramgb(blackbox) -> basis

Computes the Groebner basis of the ideal represented with the given `blackbox`.
Admissible blackboxes are subtypes of `AbstractBlackboxIdeal`.
"""
function paramgb(blackbox::T; kwargs...) where {T <: AbstractBlackboxIdeal}
    up_to_degree = get(kwargs, :up_to_degree, (Inf, Inf))
    @assert all(up_to_degree .> 0) "Total degrees must be greater than 0"
    estimate_degrees = get(kwargs, :estimate_degrees, true)
    rational_interpolator = get(kwargs, :rational_interpolator, :VanDerHoevenLecerf)
    polynomial_interpolator = get(kwargs, :polynomial_interpolator, :PrimesBenOrTiwari)
    assess_correctness = get(kwargs, :assess_correctness, false)
    ordering = get(kwargs, :ordering, Groebner.InputOrdering())
    up_to_degree_ = map(d -> isinf(d) ? div(typemax(Int), 2) : d, up_to_degree)
    ord = AbstractAlgebra.ordering(parent(blackbox))
    # If no hint for degrees is given, then try to guess degrees
    if !haskey(kwargs, :up_to_degree)
        if !estimate_degrees
            @warn "Changing the value of the parameter `estimate_degrees` from `false` to `true` since the keyword argument `up_to_degree` is not provided"
            estimate_degrees = true
        end
    end
    @info """
    Computing parametric Groebner basis up to degrees $up_to_degree
    Ordering, input / target: $ord / $(typeof(ordering))
    Rational interpolator: $rational_interpolator
    Polynomial interpolator: $polynomial_interpolator
    Estimate degrees: $estimate_degrees
    Assess correctness: $assess_correctness"""
    _paramgb(
        blackbox,
        ordering,
        up_to_degree_,
        estimate_degrees,
        assess_correctness,
        rational_interpolator,
        polynomial_interpolator
    )
end

# Progress bars
const _progressbar_color = :cyan
const _progressbar_value_color = :cyan # :light_grey
progressbar_enabled() = Logging.min_enabled_level(current_logger()) >= Logging.Info

function _paramgb(
    blackbox,
    ordering,
    up_to_degree,
    estimate_degrees,
    assess_correctness,
    rational_interpolator,
    polynomial_interpolator
)
    # Keep track of modular computation related stuff
    modular = ModularTracker(blackbox)
    # Store the state of the computation
    state = GroebnerState(blackbox, ordering)
    # Discover the shape of the groebner basis:
    # its size, and the sizes of polynomials it contains
    discover_shape!(state, modular)
    # Discover the total degrees of parameters in the coefficients. 
    if estimate_degrees
        discover_param_total_degrees!(state, modular)
    end
    # Interpolate the exponents in the parametric coefficients.
    # This uses exactly 1 prime number
    InterpolatorType = select_interpolator(rational_interpolator, polynomial_interpolator)
    @label InterpolateUsingOnePrime
    interpolate_param_exponents!(state, modular, up_to_degree, InterpolatorType)
    # Interpolate the rational coefficients of the parametric coefficients.
    # This uses the currently accumulated bases modulo different primes to
    # recover the coefficients by the means of CRT and Rational number
    # reconstruction 
    success = recover_coefficients!(state, modular, assess_correctness)
    if !success
        find_next_lucky_prime!(modular)
        @goto InterpolateUsingOnePrime
    end
    # Construct the final basis 
    basis = construct_basis(state)
    basis
end

function select_interpolator(rational_interpolator, polynomial_interpolator)
    if rational_interpolator === :VanDerHoevenLecerf
        VanDerHoevenLecerf
    else
        CuytLee
    end
end

# Discovers shape of the groebner basis of the ideal from `state`
# by specializing it at a random point (preferably, modulo a prime).
#
# If η > 0 is given, the algorithm will confirm the shape
# by additionaly specializing it at η random points.
function discover_shape!(state, modular; η=2)
    iszero(η) && (@warn "Fixing the shape of the basis from 1 point is adventurous.")
    @info "Specializing at $(1 + η) points to guess the shape of the basis.."
    blackbox = state.blackbox
    ord = state.ordering
    # Guess the shape for 1 lucky prime:
    reduce_mod_p!(blackbox, modular.ff)
    prog = ProgressUnknown(
        "# Computing specializations.. ",
        spinner=true,
        dt=0.3,
        enabled=progressbar_enabled(),
        color=_progressbar_color
    )
    @label Start
    # specialize at a random lucky point and compute GBs
    randompoints = map(_ -> randluckyspecpoint(state, modular.ff), 1:(1 + η))
    polysspecmodp = map(point -> specialize_mod_p(blackbox, point), randompoints)
    context, gb = groebner_learn(polysspecmodp[1], ordering=ord)
    state.context = context
    bases = empty(polysspecmodp)
    for i in 1:length(polysspecmodp)
        F = polysspecmodp[i]
        flag, gb = groebner_apply!(context, F)
        update!(prog, i, spinner="⌜⌝⌟⌞", valuecolor=_progressbar_value_color)
        if !flag
            @warn "Unlucky cancellation of coefficients encountered"
            @goto Start
        end
        push!(bases, gb)
    end
    # decide the "right" basis according to the majority rule
    success, basis = majorityrule(bases)
    # if not all of the bases have the same shape
    if !success
        η = 2 * (1 + η)
        @goto Start
    end
    finish!(prog)
    state.shape = basisshape(basis)
    state.basis_at_random_point[randompoints[1]] = basis
    @debug "The shape of the basis is: $(length(basis)) polynomials"
    @debug "Monomials in the basis are:" state.shape
    nothing
end

is_tdeg_interpolated_heuristic(dD, DD, dN, DN) =
    (dD + 3) < div(3DD, 4) && (dN + 3) < div(3DN, 4)

function discover_param_total_degrees!(state, modular)
    @info "Specializing at random points to guess the total degrees in parameters.."
    blackbox = state.blackbox
    ord = state.ordering
    Ru, _ = PolynomialRing(modular.ff, :u)
    K = base_ring(Ru)
    Rx = parent(blackbox)
    Ra = base_ring(Rx)
    n = length(gens(Ra))
    shift = distinct_nonzero_points(K, n)
    dilation = distinct_nonzero_points(K, n)
    shape = state.shape
    context = state.context
    # Current bounds on the total degree
    N, D = 1, 1
    npoints = N + D + 2
    # Buffers for storing evaluated and interpolated coefficients
    coeffs = Vector{Vector{Vector{elem_type(K)}}}(undef, length(shape))
    interpolated = Vector{Vector{Tuple{elem_type(Ru), elem_type(Ru)}}}(undef, length(shape))
    degrees = Vector{Vector{Tuple{Int, Int}}}(undef, length(shape))
    for i in 1:length(shape)
        interpolated[i] =
            Vector{Tuple{elem_type(Ru), elem_type(Ru)}}(undef, length(shape[i]))
        degrees[i] = Vector{Tuple{Int, Int}}(undef, length(shape[i]))
        coeffs[i] = Vector{Vector{elem_type(K)}}(undef, length(shape[i]))
        for j in 1:length(shape[i])
            coeffs[i][j] = Vector{elem_type(K)}(undef, npoints)
        end
    end
    # Interpolate coefficients on the range of J..npoints until interpolation
    # succeeds
    J = 1
    univ_points = Vector{elem_type(K)}()
    all_interpolated = false
    prog = ProgressUnknown(
        "# Computing specializations.. ",
        spinner=true,
        dt=0.3,
        enabled=progressbar_enabled(),
        color=_progressbar_color
    )
    while !all_interpolated
        all_interpolated = true
        N, D = 2N, 2D
        interpolator = CauchyInterpolator(Ru, N, D)
        npoints = N + D + 2
        @debug "Using $npoints points.."
        # TODO: points such that denominators and leading coeffs do not vanish
        univ_points = distinct_nonzero_points(K, npoints - length(univ_points), univ_points)
        x_points = map(point -> dilation .* repeat([point], n) .+ shift, univ_points)
        for i in 1:length(shape)
            for j in 1:length(shape[i])
                resize!(coeffs[i][j], npoints)
            end
        end
        for idx in J:npoints
            point = x_points[idx]
            Ip = specialize_mod_p(blackbox, point)
            flag, basis = groebner_apply!(context, Ip)
            update!(prog, idx, spinner="⌜⌝⌟⌞", valuecolor=_progressbar_value_color)
            # TODO: just select another batch of points, no need to throw
            !flag && __throw_unlucky_cancellation()
            !check_shape(shape, basis) && __throw_unlucky_cancellation()
            for i in 1:length(coeffs)
                for j in 1:length(coeffs[i])
                    coeffs[i][j][idx] = coeff(basis[i], j)
                end
            end
        end
        for i in 1:length(coeffs)
            for j in 1:length(coeffs[i])
                P, Q = interpolate!(interpolator, univ_points, coeffs[i][j])

                interpolated[i][j] = (P, Q)
                dp, dq = degree(P), degree(Q)
                degrees[i][j] = (dp, dq)
                if !is_tdeg_interpolated_heuristic(dp, N, dq, D)
                    all_interpolated = false
                end
            end
        end
        J = npoints + 1
    end
    finish!(prog)
    state.param_degrees = degrees
    @info "Success! $(npoints) points used."
    @info "The total degrees in the coefficients" state.param_degrees
    nothing
end

function interpolate_param_exponents!(
    state,
    modular,
    up_to_degree,
    ::Type{InterpolatorType}
) where {InterpolatorType}
    @info "Interpolating the exponents in parameters.."
    blackbox = state.blackbox
    ord = state.ordering
    reduce_mod_p!(blackbox, modular.ff)
    Rx = parent(blackbox)
    Ra = base_ring(Rx)
    Ru, _ = PolynomialRing(modular.ff, symbols(Ra), ordering=Nemo.ordering(Ra))
    K = base_ring(Ru)
    n = length(gens(Ra))
    shape = state.shape
    context = state.context
    # Calculate the degrees for interpolation
    Nd = up_to_degree[1] + 1
    Dd = up_to_degree[2] + 1
    total_degrees = state.param_degrees
    Nd_estimated = maximum(d -> maximum(dd -> dd[1], d), total_degrees)
    Dd_estimated = maximum(d -> maximum(dd -> dd[2], d), total_degrees)
    Nd = min(Nd, Nd_estimated)
    Dd = min(Dd, Dd_estimated)
    Nds, Dds = repeat([Nd], n), repeat([Dd], n)

    if !is_interpolation_feasible(max(Nd, Dd), K, n)
        @warn "In Prime number approach the field order might be too small" Nd Dd max(
            Nd,
            Dd
        ) * log(n) log(BigInt(order(K)))
    end

    # The current number of terms
    Nt, Dt = 1, 1
    npoints = 0
    param_exponents =
        Vector{Vector{Tuple{elem_type(Ru), elem_type(Ru)}}}(undef, length(shape))
    coeffs = Vector{Vector{Vector{elem_type(K)}}}(undef, length(shape))
    for i in 1:length(shape)
        param_exponents[i] =
            Vector{Tuple{elem_type(Ru), elem_type(Ru)}}(undef, length(shape[i]))
        coeffs[i] = Vector{Vector{elem_type(K)}}(undef, length(shape[i]))
        for j in 1:length(shape[i])
            coeffs[i][j] = Vector{elem_type(K)}(undef, npoints)
            param_exponents[i][j] = Ru(1), Ru(1)
        end
    end

    # TODO: check the feasibility of interpolation here!
    interpolator = InterpolatorType(Ru, Nd, Dd, Nds, Dds, Nt, Dt)

    J = 1
    all_interpolated = false
    @info """
    Interpolating for degrees:
    Numerator: $Nd, Denominator: $Dd"""
    prog = ProgressUnknown(
        "# Computing specializations..",
        spinner=true,
        dt=0.3,
        enabled=progressbar_enabled(),
        color=_progressbar_color
    )
    while !all_interpolated
        all_interpolated = true

        x_points = get_evaluation_points!(interpolator)
        npoints = length(x_points)
        Nt, Dt = div(npoints, 2), div(npoints, 2)
        @debug "Using $npoints points.."

        for i in 1:length(shape)
            for j in 1:length(shape[i])
                resize!(coeffs[i][j], npoints)
            end
        end

        for idx in J:npoints
            point = x_points[idx]
            Ip = specialize_mod_p(blackbox, point)
            flag, basis = groebner_apply!(context, Ip)
            update!(
                prog,
                idx,
                spinner="⌜⌝⌟⌞",
                showvalues=[(:Points, idx)],
                valuecolor=_progressbar_value_color
            )
            !flag && __throw_unlucky_cancellation()
            !check_shape(shape, basis) && __throw_unlucky_cancellation()

            for i in 1:length(shape)
                for j in 1:length(shape[i])
                    coeffs[i][j][idx] = coeff(basis[i], j)
                end
            end
        end

        for i in 1:length(shape)
            !all_interpolated && break
            for j in 1:length(shape[i])
                if total_degrees[i][j][1] > up_to_degree[1] ||
                   total_degrees[i][j][2] > up_to_degree[2]
                    continue
                end
                P, Q = interpolate!(interpolator, coeffs[i][j])
                param_exponents[i][j] = (P, Q)
                dp, dq = total_degree(P), total_degree(Q)
                tp, tq = length(P), length(Q)
                if dp < total_degrees[i][j][1] || dq < total_degrees[i][j][2]
                    all_interpolated = false
                    break
                end
                # Check that the interpolated coefficient is correct with a high
                # probability
                # c_true = coeffs_at_x0[i][j]
                # c_interpolated = evaluate(P, point_x0) // evaluate(Q, point_x0)
                # if !(c_true == c_interpolated)
                #     all_interpolated = false
                #     break
                # end
            end
        end
        J = npoints + 1
        # if npoints > 2 * (Nd + Dd) * BigInt(n) ^ max(Nd, Dd)
        #     __throw_unlucky_cancellation()
        # end
    end
    finish!(prog)
    state.param_exponents = param_exponents
    state.field_to_param_exponents[modular.ff] = state.param_exponents
    @info "Success! $(npoints) points used."
    maxDn = maximum(l -> maximum(total_degree ∘ first, l), param_exponents)
    maxDd = maximum(l -> maximum(total_degree ∘ last, l), param_exponents)
    maxTn = maximum(l -> maximum(length ∘ first, l), param_exponents)
    maxTd = maximum(l -> maximum(length ∘ last, l), param_exponents)
    @info "Basis exponents summary:
    Maximal interpolated degrees are: $maxDn for num. and $maxDd for den.
    Maximal number of interpolated terms are: $maxTn for num. and $maxTd for den.
    "
    nothing
end

function recover_coefficients!(state, modular, assess_correctness)
    @info "Recovering the coefficients.."
    reconstruct_crt!(state, modular)
    success = reconstruct_rational!(state, modular)
    if !success
        @info "Rational reconstrction failed, selecting next prime.."
        return success
    end
    if !assess_correctness
        return true
    end
    success = assess_correctness_mod_p(state, modular)
    if success
        @info "Success! Used $(length(modular.used_primes) + 1) prime in total"
    else
        @info "Correctness check failed, selecting next prime.."
    end
    success
end

function construct_basis(state)
    state.param_basis
end

function check_shape(shape, basis)
    if length(shape) != length(basis)
        return false
    end
    if map(length, shape) != map(length, basis)
        return false
    end
    # if !(shape == basisshape(basis))
    #     return false
    # end
    true
end

function majorityrule(bases::Vector{T}) where {T}
    if length(bases) == 1
        return true, first(bases)
    end
    shapes = map(basisshape, bases)
    if length(unique!(shapes)) > 1
        @warn "One of the computed bases has a different shape"
        return false, first(bases)
    end
    true, first(bases)
end

basisexponents(basis) = map(collect ∘ exponent_vectors, basis)
basisshape(basis) = map(collect ∘ monomials, basis)
basiscoeffs(basis) = map(collect ∘ coefficients, basis)
