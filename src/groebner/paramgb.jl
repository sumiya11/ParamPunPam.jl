
"""
    paramgb(polys; options...) -> basis

Computes the Groebner basis of the ideal generated by `polys`. The argument
`polys` must be an array of polynomials over a field of rational functions over
the rationals or over a finite field.

The algorithm is probabilistic and succeeds with a high probability.

## Options

Supported keyword arguments:

- `ordering`: monomial ordering. Same as for **Groebner.jl**, for example, `Lex()`
  or `DegRevLex()`.
- `up_to_degree`: Compute parameters up to a fixed total degree of numerators
  and denominators e.g., `up_to_degree=(2, 2)`. **NOTE:** If `up_to_degree` is
  specified, then a coefficient in the basis is guaranteed to be correct only if
  its total degree is smaller than or equal to `up_to_degree`. Otherwise, a
  coefficient in the basis is set to `1`.
- `rational_interpolator`: Rational function interpolation algorithm.
  Possible options are `:CuytLee` and `:VanDerHoevenLecerf` (default).
- `estimate_degrees`: If `true`, estimates the total degrees of parameters
  before starting the interpolation. Default is `true`.
- `assess_correctness`: If `true`, check that the basis is correct with high
  probability. Default is `true`. **NOTE**:

## Example

```jldoctest
using Nemo, ParamPunPam

Rparam, (a, b) = polynomial_ring(QQ, ["a", "b"])
R, (x, y, z) = polynomial_ring(fraction_field(Rparam), ["x", "y", "z"])

ParamPunPam.paramgb([a*x^2 + 1, y^2*z + (1//b)*y])
```
"""
function paramgb(polys::Vector{T}; kwargs...) where {T}
    @assert !isempty(polys) "Empty input is not supported"
    paramgb(BasicBlackboxIdeal(polys); kwargs...)
end

const _supported_kwargs = Set{Symbol}([
    :up_to_degree,
    :estimate_degrees,
    :rational_interpolator,
    :polynomial_interpolator,
    :assess_correctness,
    :ordering
])

"""
    paramgb(blackbox) -> basis

Computes the Groebner basis of the ideal represented with the given `blackbox`.
Admissible blackboxes are subtypes of `AbstractBlackboxIdeal`.
"""
function paramgb(blackbox::T; kwargs...) where {T <: AbstractBlackboxIdeal}
    for (k, _) in kwargs
        @assert k in _supported_kwargs """
        Keyword argument `$k` is not supported. 
        Supported keyword arguments are: $(join(map(string, collect(_supported_kwargs)), ", ")).

        Type `?` followed by `paramgb` to get a description of all supported
        keyword arguments."""
    end
    up_to_degree = get(kwargs, :up_to_degree, (Inf, Inf))
    @assert all(up_to_degree .> 0) "Total degrees must be greater than 0"
    estimate_degrees = get(kwargs, :estimate_degrees, true)
    rational_interpolator = get(kwargs, :rational_interpolator, :VanDerHoevenLecerf)
    @assert rational_interpolator in (:VanDerHoevenLecerf, :CuytLee)
    polynomial_interpolator = get(kwargs, :polynomial_interpolator, :PrimesBenOrTiwari)
    assess_correctness = get(kwargs, :assess_correctness, true)
    if assess_correctness && up_to_degree != (Inf, Inf)
        @debug "Turning off `assess_correctness` because `up_to_degree` was provided."
        assess_correctness = false
    end
    ordering = get(kwargs, :ordering, Groebner.InputOrdering())
    up_to_degree_ = map(d -> isinf(d) ? div(typemax(Int), 2) : d, up_to_degree)
    ord = AbstractAlgebra.internal_ordering(parent(blackbox))
    # If no hint for degrees is given, then try to guess degrees
    if !haskey(kwargs, :up_to_degree)
        if !estimate_degrees
            @warn """
            Changing the value of keyword argument `estimate_degrees` from
            `false` to `true` since the value of `up_to_degree` is not
            specified.
            """
            estimate_degrees = true
        end
    end
    _runtime_data[:npoints_degree_estimation] = 0
    _runtime_data[:npoints_interpolation] = 0
    @debug """
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
        discover_total_degrees!(state, modular, up_to_degree)
    end
    # Interpolate the exponents in the parametric coefficients.
    # This uses exactly 1 prime number
    InterpolatorType = select_interpolator(rational_interpolator, polynomial_interpolator)
    @label InterpolateUsingOnePrime
    interpolate_exponents!(state, modular, up_to_degree, InterpolatorType)
    # Interpolate the rational coefficients of the parametric coefficients. This
    # uses the currently accumulated bases modulo several different primes to
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
    # Currently, we always use PrimesBenOrTiwari for interpolating multivariate
    # polynomials
    if rational_interpolator === :VanDerHoevenLecerf
        VanDerHoevenLecerf
    else
        @assert rational_interpolator === :CuytLee
        CuytLee
    end
end

# Discovers the shape of the groebner basis of the ideal from `state` by
# specializing it at a random point (preferably, modulo a prime).
#
# If η > 0 is given, the algorithm will confirm the shape by additionaly
# specializing the basis at η random points.
function discover_shape!(state, modular; η=2)
    iszero(η) && (@warn "Discovering the shape of the basis from 1 point is adventurous.")
    @debug "Specializing at $(1 + η) points to guess the shape of the basis.."
    blackbox = state.blackbox
    ord = state.gb_ordering
    # Guess the shape for 1 lucky prime:
    reduce_mod_p!(blackbox, modular.finite_field)
    prog = ProgressUnknown(
        desc="# Computing specializations.. ",
        spinner=true,
        dt=0.3,
        enabled=is_progressbar_enabled(),
        color=_progressbar_color
    )
    @label Start
    # specialize at a random lucky point and compute GBs
    randompoints = map(_ -> randluckyspecpoint(state, modular.finite_field), 1:(1 + η))
    polysspecmodp = map(point -> specialize_mod_p(blackbox, point), randompoints)
    gb_context, gb =
        groebner_learn(polysspecmodp[1], ordering=ord, loglevel=groebner_loglevel())
    state.gb_context = gb_context
    bases = empty(polysspecmodp)
    for i in 1:length(polysspecmodp)
        F = polysspecmodp[i]
        flag, gb = groebner_apply!(gb_context, F, loglevel=groebner_loglevel())
        update!(prog, i, spinner=_progressbar_spinner, valuecolor=_progressbar_value_color)
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
    @debug "The shape of the basis is: $(length(basis)) polynomials"
    @debug "Monomials in the basis are:" state.shape
    nothing
end

const DEGREE_TOO_LARGE = -1

# Estimates the total degree of parameters of each of the coefficients of the
# Groebner basis. If the total degree of some coefficient is larger than the
# given `up_to_degree`, then the coefficient is marked with DEGREE_TOO_LARGE.
function discover_total_degrees!(state, modular, up_to_degree)
    @debug "Specializing at random points to guess the total degrees in parameters.."
    blackbox = state.blackbox
    ord = state.gb_ordering
    Ru, _ = polynomial_ring(modular.finite_field, :u)
    K = base_ring(Ru)
    Rx = parent(blackbox)
    Ra = base_ring(Rx)
    n = length(gens(Ra))
    shift = distinct_nonzero_points(K, n)
    dilation = distinct_nonzero_points(K, n)
    shape = state.shape
    gb_context = state.gb_context
    # Current bounds on the total degree
    N, D = 1, 1
    npoints = N + D + 2
    # Buffers for storing evaluated and interpolated coefficients
    coeffs = Vector{Vector{Vector{elem_type(K)}}}(undef, length(shape))
    interpolated = Vector{Vector{Tuple{elem_type(Ru), elem_type(Ru)}}}(undef, length(shape))
    degrees = Vector{Vector{Tuple{Int, Int}}}(undef, length(shape))
    @inbounds for i in 1:length(shape)
        interpolated[i] =
            Vector{Tuple{elem_type(Ru), elem_type(Ru)}}(undef, length(shape[i]))
        degrees[i] = Vector{Tuple{Int, Int}}(undef, length(shape[i]))
        coeffs[i] = Vector{Vector{elem_type(K)}}(undef, length(shape[i]))
        for j in 1:length(shape[i])
            degrees[i][j] = (DEGREE_TOO_LARGE, DEGREE_TOO_LARGE)
            coeffs[i][j] = Vector{elem_type(K)}(undef, npoints)
            interpolated[i][j] = (zero(Ru), zero(Ru))
        end
    end
    # Interpolate coefficients on the range of J..npoints until interpolation
    # succeeds
    J = 1
    univ_points = Vector{elem_type(K)}()
    all_interpolated = false
    maybe_all_interpolated = false
    prog = ProgressUnknown(
        desc="# Computing specializations.. ",
        spinner=true,
        dt=0.3,
        enabled=is_progressbar_enabled(),
        color=_progressbar_color
    )
    while !all_interpolated
        # If there is a suspision that every coefficient has been interpolated,
        # then try to use a sharper number of evaluation points
        if maybe_all_interpolated
            N, D = N + 1, D + 1
        else
            N, D = 2N, 2D
        end
        maybe_all_interpolated = true
        all_interpolated = true
        interpolator = CauchyInterpolator(Ru, N, D)
        npoints = N + D + 2
        @debug "Using $npoints points.."
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
            flag, basis = groebner_apply!(gb_context, Ip, loglevel=groebner_loglevel())
            update!(
                prog,
                idx,
                spinner=_progressbar_spinner,
                valuecolor=_progressbar_value_color
            )
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
                if maybe_all_interpolated
                    if interpolated[i][j] != (P, Q)
                        all_interpolated = false
                    end
                end
                interpolated[i][j] = (P, Q)
                dp, dq = degree(P), degree(Q)
                degrees[i][j] = (dp, dq)
                if dp >= N || dq >= D
                    maybe_all_interpolated = false
                end
            end
        end
        all_interpolated = maybe_all_interpolated && all_interpolated
        J = npoints + 1
        # The number of evaluation points must be sharp with respect to
        # `up_to_degree`
        if (up_to_degree[1] + 1) <= N && (up_to_degree[2] + 1) <= D
            break
        end
    end
    finish!(prog)
    for i in 1:length(degrees)
        for j in 1:length(degrees[i])
            dp, dq = degrees[i][j]
            if !(dp < (up_to_degree[1] + 1) && dq < (up_to_degree[1] + 1))
                degrees[i][j] = (DEGREE_TOO_LARGE, DEGREE_TOO_LARGE)
            end
        end
    end
    state.param_degrees = degrees
    _runtime_data[:npoints_degree_estimation] = npoints
    @debug "Success! $(npoints) points used."
    @debug "The total degrees in the coefficients" state.param_degrees
    nothing
end

# Interpolates the exponents of the parametric coefficients of the Groebner
# basis. Assumes that the order of the selected finite field is large enough for
# this.
function interpolate_exponents!(
    state,
    modular,
    up_to_degree,
    ::Type{InterpolatorType}
) where {InterpolatorType}
    @debug "Interpolating the exponents in parameters.."
    blackbox = state.blackbox
    ord = state.gb_ordering
    reduce_mod_p!(blackbox, modular.finite_field)
    Rx = parent(blackbox)
    Ra = base_ring(Rx)
    Ru, _ = polynomial_ring(
        modular.finite_field,
        symbols(Ra),
        internal_ordering=Nemo.internal_ordering(Ra)
    )
    K = base_ring(Ru)
    n = length(gens(Ra))
    shape = state.shape
    gb_context = state.gb_context
    # Calculate the degrees for interpolation
    Nd = up_to_degree[1]
    Dd = up_to_degree[2]
    total_degrees = state.param_degrees
    Nd_estimated = maximum(d -> maximum(dd -> dd[1], d), total_degrees)
    Dd_estimated = maximum(d -> maximum(dd -> dd[2], d), total_degrees)
    Nd = min(Nd, Nd_estimated)
    Dd = min(Dd, Dd_estimated)
    Nds, Dds = repeat([Nd], n), repeat([Dd], n)

    if !is_interpolation_feasible(max(Nd, Dd), K, n)
        @warn "In the prime number interpolation approach the field order might be too small" Nd Dd n max(
            Nd,
            Dd
        ) * log(
            n
        ) log(BigInt(order(K)))
    end

    # The current number of terms
    Nt, Dt = 1, 1
    npoints = 0
    param_exponents =
        Vector{Vector{Tuple{elem_type(Ru), elem_type(Ru)}}}(undef, length(shape))
    coeffs = Vector{Vector{Vector{elem_type(K)}}}(undef, length(shape))
    must_be_interpolated = Vector{Vector{Bool}}(undef, length(shape))
    @inbounds for i in 1:length(shape)
        param_exponents[i] =
            Vector{Tuple{elem_type(Ru), elem_type(Ru)}}(undef, length(shape[i]))
        coeffs[i] = Vector{Vector{elem_type(K)}}(undef, length(shape[i]))
        must_be_interpolated[i] = Vector{Bool}(undef, length(shape[i]))
        for j in 1:length(shape[i])
            coeffs[i][j] = Vector{elem_type(K)}(undef, npoints)
            param_exponents[i][j] = Ru(1), Ru(1)
            must_be_interpolated[i][j] = true
            if total_degrees[i][j][1] > up_to_degree[1] ||
               total_degrees[i][j][2] > up_to_degree[2]
                must_be_interpolated[i][j] = false
            end
            if total_degrees[i][j][1] == DEGREE_TOO_LARGE ||
               total_degrees[i][j][2] == DEGREE_TOO_LARGE
                must_be_interpolated[i][j] = false
            end
        end
    end

    interpolator = InterpolatorType(Ru, Nd, Dd, Nds, Dds, Nt, Dt)

    J = 1
    @debug """
    Interpolating for degrees:
    Numerator: $Nd, Denominator: $Dd"""
    prog = ProgressUnknown(
        desc="# Computing specializations..",
        spinner=true,
        dt=0.3,
        enabled=is_progressbar_enabled(),
        color=_progressbar_color
    )
    attempts = 3
    maybe_correct_basis = false
    # The outer loop is governed by Schwarz-Zippel lemma.
    # When the basis is correct with high probability, it halts.
    while !maybe_correct_basis
        # The inner loop is governed by several heuristics
        if attempts <= 0
            __throw_something_went_wrong(
                """
                Exceeded the maximum number of attempts to interpolate the basis. 
                This should not happen normally.
                Please consider submitting a Github issue."""
            )
        end
        maybe_correct_basis = true
        all_interpolated = false
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
                flag, basis = groebner_apply!(gb_context, Ip, loglevel=groebner_loglevel())
                update!(
                    prog,
                    idx,
                    spinner=_progressbar_spinner,
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
                    if !must_be_interpolated[i][j]
                        continue
                    end
                    success, P, Q = interpolate!(interpolator, coeffs[i][j])
                    if !success
                        all_interpolated = false
                        break
                    end
                    dp, dq = total_degree(P), total_degree(Q)
                    if dp < total_degrees[i][j][1] || dq < total_degrees[i][j][2]
                        all_interpolated = false
                        break
                    end
                    # Check that the interpolated coefficient is correct with a high
                    # probability
                    c_true = coeffs[i][j][end]
                    eval_point = x_points[end]
                    c_interpolated = evaluate(P, eval_point) // evaluate(Q, eval_point)
                    if !(c_true == c_interpolated)
                        all_interpolated = false
                        break
                    end
                    param_exponents[i][j] = (P, Q)
                end
            end
            J = npoints + 1
        end
        # Check that the interpolated result is correct at a random point
        random_point = distinct_nonzero_points(K, n)
        Ip = specialize_mod_p(blackbox, random_point)
        flag, gb_at_random_point =
            groebner_apply!(gb_context, Ip, loglevel=groebner_loglevel())
        !flag && __throw_unlucky_cancellation()
        @debug """
        Checking interpolated coefficients at a random points. 
        Point: $random_point
        Basis: $gb_at_random_point
        Interpolated coeffs: $param_exponents
        The number of eval. points: $npoints
        Global index: $J"""
        @inbounds for i in 1:length(shape)
            if length(gb_at_random_point[i]) != length(param_exponents[i])
                maybe_correct_basis = false
            end
            !maybe_correct_basis && break
            for j in 1:length(gb_at_random_point[i])
                if !must_be_interpolated[i][j]
                    continue
                end
                gb_coeff = coeff(gb_at_random_point[i], j)
                interpolated_num, interpolated_den = param_exponents[i][j]
                evaluated_coeff =
                    evaluate(interpolated_num, random_point) //
                    evaluate(interpolated_den, random_point)
                if gb_coeff != evaluated_coeff
                    maybe_correct_basis = false
                    break
                end
            end
        end
        attempts -= 1
    end
    finish!(prog)
    state.param_exponents = param_exponents
    state.field_to_param_exponents[modular.finite_field] = state.param_exponents
    _runtime_data[:npoints_interpolation] = npoints
    @debug "Success! $(npoints) points used."
    maxDn = maximum(l -> maximum(total_degree ∘ first, l), param_exponents)
    maxDd = maximum(l -> maximum(total_degree ∘ last, l), param_exponents)
    maxTn = maximum(l -> maximum(length ∘ first, l), param_exponents)
    maxTd = maximum(l -> maximum(length ∘ last, l), param_exponents)
    @debug """
    Basis interpolated exponents summary:
    Maximal interpolated degrees are: $maxDn for num. and $maxDd for den.
    Maximal number of interpolated terms are: $maxTn for num. and $maxTd for den.
    Points used: $npoints.
    """
    nothing
end

function recover_coefficients!(state, modular, assess_correctness)
    @debug "Recovering the coefficients.."
    reconstruct_crt!(state, modular)
    success = reconstruct_rational!(state, modular)
    if !success
        @info "Rational reconstruction failed, selecting next prime.."
        return success
    end
    if !assess_correctness
        return true
    end
    success = assess_correctness_mod_p(state, modular)
    if success
        @debug "Success! Used $(length(modular.used_primes) + 1) prime in total"
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
        @warn """
        The shape of one of the specialized Groebner bases is different from the
        others.
        """
        return false, first(bases)
    end
    true, first(bases)
end

basisexponents(basis) = map(collect ∘ exponent_vectors, basis)
basisshape(basis) = map(collect ∘ monomials, basis)
basiscoeffs(basis) = map(collect ∘ coefficients, basis)
