# Ben-or and Tiwaris via prime numbers

const _first_primes = Primes.nextprimes(2, 200)
const _prime_to_idx = Dict(
    _first_primes .=> 1:length(_first_primes)
)

mutable struct PrimesBenOrTiwari{Ring}
    # multivariate polynomial ring
    ring::Ring
    # the number of terms in the interpolant
    T::Int
    # the vector of partial degrees of the interpolant
    # Ds::Vector{Int}
    # the vector of prime numbers used in substitution
    ps::Vector{UInt}
    function PrimesBenOrTiwari(ring::Ring, T::Integer, D::Integer) where {Ring}
        @assert T >= 0
        K = base_ring(ring)
        n = nvars(ring)
        @assert (order(K) - 1) > 2T
        ps = Primes.nextprimes(2, n)
        D*log(ps[end]) >= log(UInt(order(K))) && @warn "In Prime number approach the field order might be too small" ps D order(K)
        new{Ring}(ring, T, ps)
    end
end

# 2, 3, .., pn
function startingpoint(fmbot::PrimesBenOrTiwari)
    K = base_ring(fmbot.ring)
    map(x -> K(x), fmbot.ps)
end

function factor_exponents(mi, n)
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

function interpolate!(fbot::PrimesBenOrTiwari, xs, ys)
    # check that the first degree is 0
    @assert all(isone, xs[1])
    @assert length(xs) == length(ys) == 2*fbot.T
    ω = xs[2]
    Rx = fbot.ring
    K = base_ring(Rx)
    T = fbot.T
    # construct the polynomial ys[1]z^1 + ys[2]z^2 + ... + ys[2T]z^(2T)
    # O(1)
    Rz, z = K["z"]
    sequence = Rz(ys)
    # find A/B such that A/B = sequence mod z^(2T) and degree(A) < T
    # O(M(T)logT)
    _, B, _ = Padé(sequence, z^(2T), T - 1)
    # @info "" ω B
    # @assert degree(B) == T
    # assuming this is O(T logT^k logq^m) for some k and m, 
    # where q is the order of the base field
    mi = roots(B)
    any(iszero, mi) && return zero(Rx)
    mi = map(inv, mi)
    # @assert length(mi) == T
    # @info "" mi
    # find the monomials of the interpolant,
    # O(TlogTlogq), where q is the order of the base field
    # (note that this cost covers the case where K is not a prime field)
    # (assuming ord is smooth)
    monoms = factor_exponents(mi, nvars(Rx))
    # @info "" monoms
    # find the coefficients of the interpolant
    # by solving a T×T Vandermonde system
    # O(M(T)logT).
    # t is the true number of terms
    t = min(T, length(mi))
    iszero(t) && return zero(Rx)
    coeffs = solve_transposed_vandermonde(Rz, view(mi, 1:t), view(ys, 1:t))
    Rx(coeffs, monoms)
end
