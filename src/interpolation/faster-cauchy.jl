# Cauchy interpolation using the fast polynomial gcd.

mutable struct FasterCauchy{Ring}
    ring::Ring
    N::Int
    D::Int

    function FasterCauchy(ring::Ring, N::Integer, D::Integer) where {Ring}
        @assert N >= 0 && D >= 0
        new{Ring}(ring, N, D)
    end
end

# Returns a tuple of polynomials (f, g), such that
# f(x)/g(x) = y, for all of the given xs and ys.
# O(M(n)logn), where n = max(degree(f), degree(g)).
function interpolate!(c::FasterCauchy, xs::Vector{T}, ys::Vector{T}) where {T}
    @assert length(xs) == length(ys) == c.N + c.D + 2
    R = c.ring
    z = gen(R)
    # F is the polynomial of degree at max N + D + 1
    # that passes through interpolation points,
    # O(M(n)logn)
    F = fastpolyinterpolate(R, xs, ys)
    # modulo = (z - x1)(z - x2)...(z - xn) for xi in xs,
    # O(M(n))
    modulo = producttree(z, xs)
    # r//t ≡ F (mod (z - x1)(z - x2)...(z - xn)),
    # O(M(n)logn)
    r, t, _ = Padé(F, modulo, c.N)
    #!isunit(gcd(r, t)) && throw("Cauchy interpolation fail.")
    normfactor = trailing_coefficient(t)
    divexact(r, normfactor), divexact(t, normfactor)
end
