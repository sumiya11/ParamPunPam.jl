
"""
    AbstractBlackboxIdeal

Blackbox ideals in Q(x)[y] that can be evaluated at x modulo a prime.

## Interface

Subtypes of `AbstractBlackboxIdeal` must implement the following functions:

- `length(<:AbstractBlackboxIdeal)`: the number of generators.
- `base_ring(<:AbstractBlackboxIdeal)`: original ground field.
- `parent(<:AbstractBlackboxIdeal)`: original parent ring.
- `parent_params(<:AbstractBlackboxIdeal)`: original coefficient parent ring.
- `reduce_mod_p!(<:AbstractBlackboxIdeal, p)`: reduces the generators modulo
  `p`.
- `specialize_mod_p(<:AbstractBlackboxIdeal, point)`: specializes the ideal at
  `point` and returns the generators of a specialized ideal in Zp[y].
"""
abstract type AbstractBlackboxIdeal end

@noinline __throw_unlucky_cancellation() =
    throw(AssertionError("Unlucky cancellation of coefficients!"))

mutable struct BasicBlackboxIdeal{PolyQQX} <: AbstractBlackboxIdeal
    polys::Vector{PolyQQX}
    polys_mod_p::Any

    function BasicBlackboxIdeal(polys::Vector{T}) where {T}
        @assert !isempty(polys)
        polys = filter(!iszero, polys)
        @debug "Constructing a blackbox from $(length(polys)) input polynomials"
        Rx = parent(polys[1])
        Ra = base_ring(Rx)
        if Ra isa Nemo.FracField
            Ra = base_ring(Ra)
            Rlifted, _ =
                PolynomialRing(Ra, map(string, gens(Rx)), ordering=Nemo.ordering(Rx))
            polys = liftcoeffs(polys, Rlifted)
        end
        K = base_ring(Ra)
        @assert K isa Nemo.FracField{Nemo.fmpz} "Only Nemo.QQ as a ground field is supported"
        new{eltype(polys)}(polys, nothing)
    end
end

AbstractAlgebra.base_ring(ideal::BasicBlackboxIdeal) =
    base_ring(base_ring(first(ideal.polys)))
AbstractAlgebra.parent(ideal::BasicBlackboxIdeal) = parent(first(ideal.polys))
parent_params(ideal::BasicBlackboxIdeal) = base_ring(parent(first(ideal.polys)))
Base.length(ideal::BasicBlackboxIdeal) = length(ideal.polys)

function reduce_mod_p!(ideal::BasicBlackboxIdeal, ff)
    @debug "Reducing modulo $(ff).."
    ideal.polys_mod_p = map(
        poly -> map_coefficients(f -> map_coefficients(c -> ff(c), f), poly),
        ideal.polys
    )
    nothing
end

function specialize_mod_p(ideal::BasicBlackboxIdeal, point)
    polys_mod_p = ideal.polys_mod_p
    for poly in polys_mod_p
        if iszero(evaluate(leading_coefficient(poly), point))
            __throw_unlucky_cancellation()
        end
    end
    map(f -> map_coefficients(c -> evaluate(c, point), f), polys_mod_p)
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
