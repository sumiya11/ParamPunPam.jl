
mutable struct GroebnerState{Blackbox, FF, PolyFF, PolyFracQQ, OrderingGb}
    # original polynomials over Q(a)
    blackbox::Blackbox
    shape::Vector{Vector{PolyFF}}
    # total degrees of the parameters of the Groebner basis
    param_degrees::Vector{Vector{Tuple{Int, Int}}}
    # exponents of the parameters of the Groebner basis
    param_exponents::Vector{Vector{Tuple{PolyFF, PolyFF}}}
    param_coeffs_crt::Vector{Vector{Tuple{Vector{BigInt}, Vector{BigInt}}}}
    field_to_param_exponents::Dict{FF, Vector{Vector{Tuple{PolyFF, PolyFF}}}}
    # fully reconstructed basis
    param_basis::Vector{PolyFracQQ}
    # GB computation helpers
    gb_context::Any
    gb_ordering::OrderingGb

    function GroebnerState(blackbox::Blackbox, ord::Ord) where {Blackbox <: AbstractBlackboxIdeal, Ord}
        Rx = parent(blackbox)
        Ra = parent_params(blackbox)
        params = gens(Ra)
        polyvars = gens(Rx)
        K = base_ring(Ra)
        @debug "Given $(length(blackbox)) functions in $K($(join(repr.(params),", ")))[$(join(repr.(polyvars),", "))]"
        PolyFF = Nemo.fpMPolyRingElem
        PolyFracQQ = AbstractAlgebra.Generic.MPoly{AbstractAlgebra.Generic.FracFieldElem{Nemo.QQMPolyRingElem}}
        FF = Nemo.fpField
        new{Blackbox, FF, PolyFF, PolyFracQQ, typeof(ord)}(
            blackbox,
            Vector{Vector{Nemo.fpMPolyRingElem}}(),
            Vector{Vector{Tuple{Int, Int}}}(),
            Vector{Vector{Tuple{PolyFF, PolyFF}}}(),
            Vector{Vector{Tuple{Vector{BigInt}, Vector{BigInt}}}}(),
            Dict{FF, Vector{Vector{Tuple{PolyFF, PolyFF}}}}(),
            Vector{PolyFracQQ}(),
            nothing,
            ord
        )
    end
end

function randluckyspecpoint(state::GroebnerState, K)
    # TODO: this is not correct!
    Ra = parent_params(state.blackbox)
    [rand(K) for _ in 1:nvars(Ra)]
end

function assess_correctness_mod_p(state, modular)
    # NOTE: in this state, ideally, this should always pass, since the modulo
    # we use to check correctness is the same as the one used to compute the
    # basis
    ff = modular.ff
    reduce_mod_p!(state.blackbox, ff)
    point = randluckyspecpoint(state, ff)
    @debug "Checking correctness at $point in $ff"
    generators_zp = specialize_mod_p(state.blackbox, point)
    R_zp = parent(first(generators_zp))
    basis_specialized_coeffs = map(
        f -> map(
            c ->
                evaluate(map_coefficients(cc -> ff(cc), numerator(c)), point) //
                evaluate(map_coefficients(cc -> ff(cc), denominator(c)), point),
            collect(coefficients(f))
        ),
        state.param_basis
    )
    param_basis_specialized = map(
        i -> R_zp(basis_specialized_coeffs[i], collect(exponent_vectors(state.param_basis[i]))),
        1:length(basis_specialized_coeffs)
    )
    @debug "Evaluated basis" param_basis_specialized
    @debug "Evaluated generators" generators_zp
    if !isgroebner(param_basis_specialized)
        @debug "Not a basis!"
        return false
    end
    inclusion = normalform(param_basis_specialized, generators_zp)
    @debug "Inclusion in correctness assessment" inclusion
    all(iszero, inclusion)
end
