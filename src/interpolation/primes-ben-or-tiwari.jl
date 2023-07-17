# Ben-or and Tiwari via prime numbers

# If n, D are the number of variables and the total degree, respectively, and K
# is the order of the ground field, then the interpolation succeeds when 
#   p_n^D < K
# or, somewhat equivalently,
#   D log n < log K

const _first_primes = Primes.nextprimes(2, 200)
const _prime_to_idx = Dict(
    _first_primes .=> 1:length(_first_primes)
)

is_interpolation_feasible(D, K, n) = D*log(_first_primes[n]) < log(BigInt(order(K)))

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

function factor_exponents(mi::Vector{T}, n::Integer) where {T}
    exps = map(m -> Primes.factor(Dict, UInt(data(m))), mi)
    monoms = Vector{Vector{UInt}}(undef, length(exps))
    for (i, e) in enumerate(exps)
        monoms[i] = zeros(UInt, n)
        for (p, e) in e
            if !haskey(_prime_to_idx, p)
                for j in 1:length(exps)
                    monoms[j] = zeros(UInt, n)
                end
                return monoms
            end
            idx = _prime_to_idx[p]
            if idx > n
                for j in 1:length(exps)
                    monoms[j] = zeros(UInt, n)
                end
                return monoms
            end
            monoms[i][idx] = e
        end
    end
    monoms
end

function interpolate!(bot::PrimesBenOrTiwari, xs, ys)
    # check that the first degree is 0
    @assert all(isone, xs[1])
    @assert length(xs) == length(ys) == 2*bot.T
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
    mi = roots(B)
    any(iszero, mi) && return one(Rx)
    mi = map(inv, mi)
    # find the monomials of the interpolant,
    # O(TlogTlogq), where q is the order of the base field
    # (note that this cost covers the case where K is not a prime field)
    # (assuming ord is smooth)
    monoms = factor_exponents(mi, nvars(Rx))
    # find the coefficients of the interpolant
    # by solving a T×T Vandermonde system
    # O(M(T)logT).
    # t is the true number of terms
    t = min(T, length(mi))
    iszero(t) && return zero(Rx)
    coeffs = solve_transposed_vandermonde(Rz, view(mi, 1:t), view(ys, 1:t))
    Rx(coeffs, monoms)
end
