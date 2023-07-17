
mutable struct GroebnerState{BB}#, D, E, C}
    # original polynomials over Q(a)
    # (`a` are transcendental parameters)
    blackbox::BB
    shape::Any#::Vector{Vector{S}}
    coeffs_at_random_point::Dict{Any, Any}
    # total degrees of the parameters `a` of the Groebner basis
    param_degrees::Any#::Vector{Vector{D}}
    # exponents of the parameters `a` of the Groebner basis
    param_exponents::Any#::Vector{Vector{E}}
    # coefficients in the parameters `a` of the Groebner basis
    # (in Q)
    param_coeffs::Any#::Vector{Vector{C}}
    field_to_polys::Any
    graph::Any

    function GroebnerState(blackbox::T) where {T <: AbstractBlackboxIdeal}
        Rx = parent(blackbox)
        Ra = parent_params(blackbox)
        params = gens(Ra)
        polyvars = gens(Rx)
        K = base_ring(Ra)
        @logmsg LogLevel(0) "Given $(length(blackbox)) functions in $K($(join(repr.(params),", ")))[$(join(repr.(polyvars),", "))]"
        new{T}(blackbox, nothing, Dict(), nothing, nothing, nothing, Dict(), nothing)
    end
end

function randluckyspecpoint(state::GroebnerState, K)
    # TODO: this is not correct!
    Ra = parent_params(state.blackbox)
    [rand(K) for _ in 1:nvars(Ra)]
end
