# Product of all (z - x) for x in xs[i..j] 
function _producttree(z, xs, i, j)
    i == j && return z - xs[i]
    m = div(i + j, 2)
    _producttree(z, xs, i, m) * _producttree(z, xs, m + 1, j)
end

# computes the product of all (z - x) for x in xs
# using the tree product algorithm.
# O(M(n)), n = length(xs)
function producttree(z, xs)
    _producttree(z, xs, 1, length(xs))
end

treeroot(tree) = tree[end][1]
treedepth(tree) = length(tree)
treebase(tree) = length(first(tree))

# Build the tree of products of (z - x) for x in xs and returns it.
# Leaves of the tree are (z - x) and the root is the product,
# O(M(n)), n = length(xs)
function buildproducttree(z::T, xs) where {T}
    n = length(xs)
    @assert n > 0  
    npow = nextpow(2, n)
    k = round(Int, log(2, npow)) + 1
    tree = Vector{Vector{T}}(undef, k)
    tree[1] = Vector{T}(undef, npow)
    @inbounds for i in 1:n
        tree[1][i] = z - xs[i]
    end
    @inbounds for i in n+1:npow
        tree[1][i] = one(z)
    end
    @inbounds for i in 2:k
        nel = 2^(k-i)
        tree[i] = Vector{T}(undef, nel)
        for j in 1:nel
            tree[i][j] = tree[i-1][2*j-1] * tree[i-1][2*j]
        end
    end
    tree
end

function _remindertree!(f, rtree, ptree, depth, idx)
    iszero(f) && return rtree
    if depth == 0
        rtree[idx] = coeff(mod(f, ptree[1][idx]), 0)
        return rtree
    end
    l, r = 2*idx - 1, 2*idx
    @inbounds r0 = mod(f, ptree[depth][l])
    @inbounds r1 = mod(f, ptree[depth][r])
    _remindertree!(r0, rtree, ptree, depth - 1, l)
    _remindertree!(r1, rtree, ptree, depth - 1, r)
    rtree
end

# Given a polynomial f and a product tree ptree
# built from (z - x) for x in xs,
# returns the vector of reminders (f mod (z - x)) for x in xs.
# O(M(n)logn), if n = degree(f) = O(length(xs))
function remindertree(f, ptree)
    K = base_ring(parent(f))
    rtree = zeros(K, treebase(ptree))
    _remindertree!(f, rtree, ptree, treedepth(ptree) - 1, 1)
end

# Solves the system
# | 1       1     ...     1    | | x1 |   | a1 |
# | v1      v2    ...     vT   | | x2 |   | a2 |
# | .       .     ...     .    | | .  | = | .  |
# | vn^k    vn^k  ...     vn^k | | xn |   | an |
# where k is n-1
# O(M(n)logn)
# Notations are taken from 
#   "Improved Sparse Multivariate Polynomial Interpolation Algorithms",
#   Erich Kaltofen and Lakshman Yagati
function solve_transposed_vandermonde(Rz, vi, ai)
    @assert length(vi) == length(ai)
    n = length(vi)
    z = gen(Rz)
    # Dz = a1^z^n + a2^z^n-1 + ... + an^z
    # O(1)
    Dz = Nemo.shift_left(Rz(reverse(ai)), 1)
    # Bz = (z - v1)(z - v2)...(z - vn)
    # O(M(n))

    ptree = buildproducttree(z, vi)
    Bz = treeroot(ptree)
    ∂B = derivative(Bz)
    # αi = ∏(vi - vj) for all j ≠ i
    # O(M(n)logn)
    αi = remindertree(∂B, ptree)
    # O(M(n))
    Qn = Nemo.shift_right(Bz * Dz, n + 1)
    # Qnvi = Qn evaluated at vi for all i
    # O(M(n)logn)
    Qnvi = remindertree(Qn, ptree)
    # O(n)
    xi = map(i -> Qnvi[i] * inv(αi[i]), 1:n)
    xi
end

function _lagrangetree(z, ys, ptree, depth, idx)
    R = parent(z)
    if depth == 0
        return R(ys[idx])
    end
    l, r = 2*idx - 1, 2*idx
    r0 = _lagrangetree(z, ys, ptree, depth - 1, l)
    r1 = _lagrangetree(z, ys, ptree, depth - 1, r)
    r0 * ptree[depth][r] + r1 * ptree[depth][l]
end

function lagrangetree(z, ys, ptree)
    _lagrangetree(z, ys, ptree, treedepth(ptree) - 1, 1) 
end

# Returns a unique univariate polynomial f in the ring R,
# such that f(x) = y for all (x, y) in xs, ys;
# O(M(n)logn), where n = O(length(xs))
function fastpolyinterpolate(R, xs, ys)
    @assert length(xs) == length(ys)
    z = gen(R)
    # O(M(n))
    ptree = buildproducttree(z, xs)
    m = treeroot(ptree)
    dm = derivative(m)
    # O(M(n)logn)
    si = remindertree(dm, ptree)
    ysi = zeros(base_ring(R), nextpow(2, length(ys)))
    for i in 1:length(ys)
        ysi[i] = ys[i]*inv(si[i])
    end
    # O(M(n)logn)
    lagrangetree(z, ysi, ptree)
end

# Adapted from
#   "Solving Systems of Non-Linear Polynomial Equations Faster"
#   by Canny, Kaltofen, and Yagati, 1989
#
# Currently, another thing is implemented
#
# Given a multivariate polynomial
#   f(x1..xn) = a1m1 + a2m2 + ... + atmt 
# and a list of evaluation points
#   ωs = [(1,...,1),..,(ω1^i,...,ωn^i),..,(ω1^T,..., ωn^{T-1})]
# Computes and returns bi = f(ωi) for ωi in ωs in O(M(T)log t) 
function fast_multivariate_evaluate(f, ωs)
    #     |1 v1 ... v1^{T-1}|
    # V = |1 v2 ... v2^{T-1}|
    #     |     ...
    #     |1 vt ... vt^{T-1}|
    # a = (a1...at), b = (b0...bT-1)
    # Then
    #   b = V^T a  (1)
    # Let a' = V \ a (solve this in M(t)log t)
    # Then (1) becomes 
    #   b = (V^T V) a'
    # The product of Hankel matrix (V^T V) and a'
    # can be computed in O(M(t)), and the entries of
    # V^T V can be computed in O(M(T)log t)
    # a = collect(coefficients(f))
    map(ω -> evaluate(f, ω), ωs)
end
