# Ben-or and Tiwari via Kronecker substitution

# If n, D are the number of variables and the total degree, respectively, and K
# is the order of the ground field, then the interpolation succeeds when 
#   D^n < K

"""
    KronBenOrTiwari

An object for interpolating multivariate polynomials.
Uses Ben-or and Tiwari algorithm with the Kronecker substitution approach.

## Usage example

```julia
using Nemo
R, (x1, x2, x3) = Nemo.Native.GF(2^62 + 135)["x1", "x2", "x3"]

poly = x1 * x2 + x2 * x3

# the number of terms, total degree
T, D = 2, 2

interpolator = ParamPunPam.KronBenOrTiwari(R, T, D)
ω = ParamPunPam.startingpoint(interpolator)

# evaluations
xs = map(i -> ω .^ i, 0:(2T - 1))
ys = map(x -> evaluate(poly, x), xs)

success, interpolated = ParamPunPam.interpolate!(interpolator, xs, ys)

@assert success && interpolated == poly
```
"""
mutable struct KronBenOrTiwari{Ring}
    # multivariate polynomial ring
    ring::Ring
    # the number of terms in the interpolant
    T::Int
    # the vector of partial degrees of the interpolant
    Ds::Vector{Int}
    # the vector of degrees used in Kronecker substitution
    Dsubs::Vector{BigInt}
    function KronBenOrTiwari(ring, T, D)
        KronBenOrTiwari(ring, T, repeat([D], nvars(ring)))
    end
    function KronBenOrTiwari(ring::Ring, T::Integer, Ds::Vector{<:Integer}) where {Ring}
        @assert T >= 0
        @assert all(>=(0), Ds)
        @assert length(Ds) == nvars(ring)
        Dsubs = subsdegrees(Ds)
        K = base_ring(ring)
        @assert (order(K) - 1) > 2T
        if Dsubs[end] >= order(K)
            @warn "In Kronecker approach the field order might be too small" Dsubs order(K)
        end
        new{Ring}(ring, T, Ds, Dsubs)
    end
end

# Returns [1, D1, D1D2, ..., D1...Dn] to be used in Kronecker substitution.
# (probably, adding 1 to each entry)
function subsdegrees(Ds::Vector{<:Integer})
    ans = Vector{BigInt}(undef, length(Ds) + 1)
    ans[1] = one(BigInt)
    for i in 1:length(Ds)
        ans[i + 1] = ans[i] * (Ds[i] + 1)
    end
    ans
end

subsbackward(bot::KronBenOrTiwari, monoms::Vector{I}) where {I} = map(m -> subsbackward(bot, m), monoms)
function subsbackward(bot::KronBenOrTiwari, monom::I) where {I}
    Dsubs = bot.Dsubs
    n = length(Dsubs) - 1
    ans = Vector{Int}(undef, n)
    for i in n:-1:1
        ans[i] = div(monom, Dsubs[i])
        monom = monom - ans[i] * Dsubs[i]
    end
    ans
end

function startingpoint(bot::KronBenOrTiwari)
    # Do Kronecker substitution
    Dsubs = bot.Dsubs
    K = base_ring(bot.ring)
    g = randomgenerator(K)
    map(i -> g^Dsubs[i], 1:(length(Dsubs) - 1))
end

function interpolate!(bot::KronBenOrTiwari, xs, ys)
    # check that the first degree is 0
    @assert all(isone, xs[1])
    @assert length(xs) == length(ys) == 2 * bot.T
    ω = xs[2]
    Rx = bot.ring
    K = base_ring(Rx)
    T = bot.T
    # construct the polynomial ys[1]z^1 + ys[2]z^2 + ... + ys[2T]z^(2T)
    # O(1)
    Rz, z = K["z"]
    sequence = Rz(ys)
    # find A/B such that A/B = sequence mod z^(2T) and degree(A) < T
    # O(M(T)logT)
    _, B, _ = Padé(sequence, z^(2T), T - 1)
    # assuming this is O(T logT^k logq^m) for some k and m, 
    # where q is the order of the base field
    mi = map(inv, roots(B))
    # find the monomials of the interpolant,
    # O(TlogTlogq), where q is the order of the base field
    # (note that this cost covers the case where K is not a prime field)
    # (assuming ord is smooth)
    buf = DiscreteLogBuffers(PrecomputedField(K))
    monoms = map(m -> discrete_log(ω[1], m, buf), mi)
    # find the coefficients of the interpolant
    # by solving a T×T Vandermonde system
    # O(M(T)logT).
    # t is the true number of terms
    t = min(T, length(mi))
    if iszero(t)
        return true, Rx(0)
    end
    coeffs = solve_transposed_vandermonde(Rz, view(mi, 1:t), view(ys, 1:t))
    true, Rx(coeffs, subsbackward(bot, monoms))
end
