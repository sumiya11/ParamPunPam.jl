# Ben-or and Tiwari via prime numbers

# If n, D are the number of variables and the total degree, respectively, and K
# is the order of the ground field, then the interpolation succeeds when 
#   p_n^D < K
# or, somewhat equivalently,
#   D log n < log K

const _first_primes = Primes.nextprimes(2, 200)
const _prime_to_idx = Dict(_first_primes .=> 1:length(_first_primes))

is_interpolation_feasible(D, K, n) = D * log(_first_primes[n]) < log(BigInt(order(K)))

"""
    PrimesBenOrTiwari

An object for interpolating multivariate polynomials.
Uses Ben-or and Tiwari algorithm with the prime-number approach.

## Usage example

```julia
using Nemo
R, (x1, x2, x3) = Nemo.Native.GF(2^62 + 135)["x1", "x2", "x3"]

poly = x1 * x2 + x2 * x3

# the number of terms, total degree
T, D = 2, 2

interpolator = ParamPunPam.PrimesBenOrTiwari(R, T, D)
ω = ParamPunPam.startingpoint(interpolator)

# evaluations
xs = map(i -> ω .^ i, 0:(2T - 1))
ys = map(x -> evaluate(poly, x), xs)

success, interpolated = ParamPunPam.interpolate!(interpolator, xs, ys)

@assert success && interpolated == poly
```
"""
mutable struct PrimesBenOrTiwari{Ring}
    # multivariate polynomial ring
    ring::Ring
    # the number of terms in the interpolant
    T::Int
    # the vector of prime numbers used in substitution
    ps::Vector{UInt}
    function PrimesBenOrTiwari(ring::Ring, T::Integer, D::Integer) where {Ring}
        @assert T >= 0
        K = base_ring(ring)
        n = nvars(ring)
        @assert (order(K) - 1) > 2T "The number of terms is too large: cannot interpolate"
        ps = _first_primes[1:n]
        new{Ring}(ring, T, ps)
    end
end

# Returns a point [2, 3, .., pn]
function startingpoint(bot::PrimesBenOrTiwari)
    K = base_ring(bot.ring)
    map(x -> K(x), bot.ps)
end

function factor_with_known_factors(e::T, factors::Vector{T}) where {T}
    exps = zeros(UInt, length(factors))
    @inbounds for i in 1:length(factors)
        r = zero(T)
        e_next = e
        j = 0
        while iszero(r)
            e = e_next
            e_next, r = divrem(e, factors[i])
            j += 1
        end
        exps[i] = j - 1
    end
    isone(e), exps
end

function factor_exponents(mi::Vector{T}, n::Integer, factors::Vector{U}) where {T, U}
    exps = Vector{Vector{UInt}}(undef, length(mi))
    i = 1
    while i <= length(mi)
        flag, factorization = factor_with_known_factors(UInt(data(mi[i])), factors)
        if !flag
            break
        end
        exps[i] = factorization
        i += 1
    end
    if i != length(mi) + 1
        resize!(exps, i - 1)
    end
    exps
end

function interpolate!(bot::PrimesBenOrTiwari, xs, ys)
    # check that the first degree is 0
    @assert all(isone, xs[1])
    @assert length(xs) == length(ys) == 2 * bot.T
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
    mi = Nemo.roots(B)
    any(iszero, mi) && return false, one(Rx)
    mi = map(inv, mi)
    # find the monomials of the interpolant,
    # O(TlogTlogq), where q is the order of the base field
    # (note that this cost covers the case where K is not a prime field)
    # (assuming ord is smooth)
    monoms = factor_exponents(mi, nvars(Rx), bot.ps)
    # find the coefficients of the interpolant
    # by solving a T×T Vandermonde system
    # O(M(T)logT).
    # t is the true number of terms
    t = min(T, length(monoms))
    success = length(monoms) == length(mi)
    success == success && (!iszero(t) || length(monoms) == length(mi))
    (!success || iszero(t)) && return success, zero(Rx)
    coeffs = solve_transposed_vandermonde(Rz, view(mi, 1:t), view(ys, 1:t))
    interpolated = Rx(coeffs, monoms)
    success, interpolated
end
