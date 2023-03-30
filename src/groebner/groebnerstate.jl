
mutable struct GroebnerState{T, F}#, D, E, C}
    # original polynomials over Q(a)
    # (`a` are transcendental parameters)
    polys::Vector{T}
    polys_fracfree::Vector{F}
    shape#::Vector{Vector{S}}
    # total degrees of the parameters `a` of the Groebner basis
    param_degrees#::Vector{Vector{D}}
    # exponents of the parameters `a` of the Groebner basis
    param_exponents#::Vector{Vector{E}}
    # coefficients in the parameters `a` of the Groebner basis
    # (in Q)
    param_coeffs#::Vector{Vector{C}}

    function GroebnerState(polys::Vector{T}) where {T}
        Rx = parent(first(polys))
        Ra = base_ring(base_ring(Rx))
        params = gens(Ra)
        polyvars = gens(Rx)
        K = base_ring(Ra)
        @info "Given $(length(polys)) functions in K($(join(repr.(params),", ")))[$(join(repr.(polyvars),", "))]"
        # Remove denominators from the input by lifting it to a polynomial ring
        Rlifted, _ = PolynomialRing(Ra, map(string, gens(Rx)))
        polys_fracfree = liftcoeffs(polys, Rlifted)
        # @info "Lifting to polynomials in K[$params][$polyvars].."
        new{T, eltype(eltype(polys_fracfree))}(
            polys, polys_fracfree, 
            nothing, nothing, nothing, nothing
        )
    end
end

function liftcoeffs(polys, newring)
    cfs = map(collect ∘ coefficients, polys)
    cfs = map(c -> map(numerator, c .* lcm(map(denominator, c))), cfs)
    newpolys = Vector{elem_type(newring)}(undef, length(polys))
    for i in 1:length(newpolys)
        G = lcm(map(denominator, collect(coefficients(polys[i]))))
        newpolys[i] = map_coefficients(c -> numerator(c * G), polys[i])
    end
    newpolys
end

function randluckyspecpoint(state::GroebnerState, K)
    Rx = parent(first(state.polys_fracfree))
    Ra = base_ring(Rx)
    [rand(K) for _ in 1:nvars(Ra)]
end
