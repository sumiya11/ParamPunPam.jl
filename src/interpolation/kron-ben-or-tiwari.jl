# Ben-or and Tiwaris via Kronecker substitution

mutable struct KronBenOrTiwari{Ring}
    # multivariate polynomial ring
    ring::Ring
    # the number of terms in the interpolant
    T::Int
    # the vector of partial degrees of the interpolant
    Ds::Vector{Int}
    # the vector of degrees used in Kronecker substitution
    Dsubs::Vector{BigInt}
    function KronBenOrTiwari(ring::Ring, T::Integer, Ds::Vector{<:Integer}) where {Ring}
        @assert T >= 0
        @assert all(>=(0), Ds)
        @assert length(Ds) == nvars(ring)
        Dsubs = subsdegrees(Ds)
        K = base_ring(ring)
        @assert (order(K) - 1) > 2T
        Dsubs[end] >= order(K) && @warn "In Kronecker substitution the field order might be too small" Dsubs order(K)
        new{Ring}(ring, T, Ds, Dsubs)
    end
end

# Returns [1, D1, D1D2, ..., D1...Dn] to be used
# in Kronecker substitution.
# (probably, adding 1 to each entry)
function subsdegrees(Ds::Vector{<:Integer})
    ans = Vector{BigInt}(undef, length(Ds) + 1)
    ans[1] = one(BigInt)
    for i in 1:length(Ds)
        ans[i+1] = ans[i]*(Ds[i] + 1)
    end
    ans
end

subsbackward(fmbot::KronBenOrTiwari, monoms::Vector{I}) where {I} = map(m -> subsbackward(fmbot, m), monoms)
function subsbackward(fmbot::KronBenOrTiwari, monom::I) where {I}
    Dsubs = fmbot.Dsubs
    n = length(Dsubs) - 1
    ans = Vector{Int}(undef, n)
    for i in n:-1:1
        ans[i] = div(monom, Dsubs[i])
        monom = monom - ans[i]*Dsubs[i]
    end
    ans
end

function startingpoint(fmbot::KronBenOrTiwari)
    # Do Kronecker substitution
    Dsubs = fmbot.Dsubs
    K = base_ring(fmbot.ring)
    g = randomgenerator(K)
    map(i -> g^Dsubs[i], 1:length(Dsubs)-1)
end

function interpolate!(fbot::KronBenOrTiwari, blackbox)
    T = fbot.T
    # the base of the geometric sequence ω^1, ω^2, ...
    ω = startingpoint(fbot)
    # generate the sequence
    # O(TlogT)
    ωs = map(i -> ω .^ i, 0:2T-1)
    # evaluate the blackbox function at the sequence
    # O(LT), where L is the cost of evaluating the blackbox function
    ys = map(blackbox, ωs)
    interpolate!(fbot, ωs, ys)
end

function interpolate!(fbot::KronBenOrTiwari, xs, ys)
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
    mi = map(inv, roots(B))
    # @assert length(mi) == T
    # @info "" mi
    # find the monomials of the interpolant,
    # O(TlogTlogq), where q is the order of the base field
    # (note that this cost covers the case where K is not a prime field)
    # (assuming ord is smooth)
    buf = DiscreteLogBuffers(PrecomputedField(K))
    monoms = map(m -> discrete_log(ω[1], m, buf), mi)
    # @info "" monoms
    # find the coefficients of the interpolant
    # by solving a T×T Vandermonde system
    # O(M(T)logT).
    # t is the true number of terms
    t = min(T, length(mi))
    coeffs = solve_transposed_vandermonde(Rz, view(mi, 1:t), view(ys, 1:t))
    Rx(coeffs, subsbackward(fbot, monoms))
end
