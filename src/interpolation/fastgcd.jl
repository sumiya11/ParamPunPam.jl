# Fast polynomial gcd algorithm.
# The notation and step numbering are taken from
# Algorithm 11.4, Modern Computer Algebra, by Gathen and Gerhard

# Given 2×2 matrix A and vector x, returns the matrix-vector product Ax.
# Input is not modified and output is not shared.
function matvec2by1(A::Tuple{V, V}, x::V) where {V<:Tuple{T, T}} where {T}
    R = parent(x[1])
    ans = (R(), R())
    tmp = R()
    ##
    Nemo.mul!(ans[1], A[1][1], x[1])
    Nemo.mul!(tmp, A[1][2], x[2])
    Nemo.add!(ans[1], ans[1], tmp)
    ##
    Nemo.mul!(ans[2], A[2][1], x[1])
    Nemo.mul!(tmp, A[2][2], x[2])
    Nemo.add!(ans[2], ans[2], tmp)
    ##
    ans
end

# Given two matrices A and B, returns the matrix product AB.
# Input is not modified and output is not shared.
function matmul2by2(A::Tup, B::Tup) where {Tup<:Tuple{Tuple{T, T}, Tuple{T, T}}} where {T}
    R = parent(A[1][1])
    ans = ((R(), R()), (R(), R()))
    tmp = R()
    ##
    Nemo.mul!(ans[1][1], A[1][1], B[1][1])
    Nemo.mul!(tmp, A[1][2], B[2][1])
    Nemo.add!(ans[1][1], ans[1][1], tmp)
    ##
    Nemo.mul!(ans[1][2], A[1][1], B[1][2])
    Nemo.mul!(tmp, A[1][2], B[2][2])
    Nemo.add!(ans[1][2], ans[1][2], tmp)
    ##
    Nemo.mul!(ans[2][1], A[2][1], B[1][1])
    Nemo.mul!(tmp, A[2][2], B[2][1])
    Nemo.add!(ans[2][1], ans[2][1], tmp)
    ##
    Nemo.mul!(ans[2][2], A[2][1], B[1][2])
    Nemo.mul!(tmp, A[2][2], B[2][2])
    Nemo.add!(ans[2][2], ans[2][2], tmp)
    ##
    ans
end

# Returns f ⥣ k, given as f quo x^(n - k) 
# (here, n = degree(f))
function ⥣(f, k)
    (k < 0) && return zero(f)
    if degree(f) < k
        Nemo.shift_left(f, k - degree(f))
    else
        Nemo.shift_right(f, degree(f) - k)
    end
end

_gcd_basecase_threshold() = 2^4

function _direct_eea(g, f, k)
    @assert degree(g) >= degree(f)
    R = parent(g)
    V = (one(R), zero(R), g)
    U = (zero(R), one(R), f)
    a, b, c = R(), R(), R()
    # in Nemo, degree(0) is -1
    while degree(U[3]) > k
        q = div(V[3], U[3])
        #
        Nemo.mul!(a, q, U[1])
        Nemo.mul!(b, q, U[2])
        Nemo.mul!(c, q, U[3])
        #
        T = V .- (a, b, c)
        V = U
        U = T
    end
    U[3], ((V[1], V[2]), (U[1], U[2]))
end

# Fast polynomial gcd.
# Assumes the degree of r0 is greater than the degree of r1
function _fastgcd(r0, r1, k)
    @assert degree(r0) > degree(r1)
    n0, n1 = degree(r0), degree(r1)
    if iszero(r1) || k < n0 - n1
        return zero(r0), ((one(r0), zero(r0)), (zero(r0), one(r0)))
    end
    if n0 < _gcd_basecase_threshold()
        return _direct_eea(r0, r1, -(k - n0 + 1))
    end
    # first recursive call
    d = div(k, 2)
    jm1, R = _fastgcd(r0 ⥣ 2d, r1 ⥣ (2d - (n0 - n1)), d)
    rjm1, rj = matvec2by1(R, (r0, r1))
    _, nj = degree(rjm1), degree(rj)
    if iszero(rj) || k < n0 - nj
        return jm1, R
    end
    qj = div(rjm1, rj)
    rjp1 = rjm1 - qj*rj
    rhojp1 = leading_coefficient(rjp1)
    # in Nemo, the leading_coefficient of 0 is 0;
    # we want it to be 1 here.
    iszero(rhojp1) && (rhojp1 = one(rhojp1))
    rjp1 = divexact(rjp1, rhojp1)
    njp1 = degree(rjp1)
    # second recursive call
    d⁺ = k - (n0 - nj)
    hmj, S = _fastgcd(rj ⥣ 2d⁺, rjp1 ⥣ (2d⁺ - (nj - njp1)), d⁺)
    Qj = ((zero(r0), one(r0)), (parent(qj)(inv(rhojp1)), -qj*inv(rhojp1)))
    hmj + (jm1 + 1), matmul2by2(S, matmul2by2(Qj, R))
end

function standardize(g, f)
    !ismonic(g) && (g = divexact(g, leading_coefficient(g)))
    !ismonic(f) && (f = divexact(f, leading_coefficient(f)))
    if degree(g) < degree(f)
        g, f = f, g
    end
    if degree(g) == degree(f)
        g, f = f, g - f 
    end
    @assert degree(g) > degree(f)
    g, f
end

function fastgcd(g, f)
    (iszero(g) || iszero(f)) && (return one(g))
    g, f = standardize(g, f)
    _, R = _fastgcd(g, f, degree(g))
    h = first(matvec2by1(R, (g, f)))
    divexact(h, leading_coefficient(h))
end

# given (polynomials) g and f (|f| >= |g|),
# computes and returns a single row from the EEA algorithm (r, t, s), 
# such that r = t*g + s*f, |r| < k, where |r| is the maximal possible
function fastconstrainedEEA(g, f, k)
    @assert degree(g) > degree(f)
    _, R = _fastgcd(g, f, degree(g) - k - 1)
    ri, rj = matvec2by1(R, (g, f))
    if degree(ri) <= k
        t, s = R[1]
        return ri, t, s 
    else
        t, s = R[2]
        return rj, t, s
    end
end

# Given (polynomials) f and g (|f| <= |g|),
# computes and returns a single row from the EEA algorithm (r, t, s), 
# such that r = t*f + s*g, |r| < k, where |r| is maximal possible.
# O(M(T)logT), where T = max(degree(f), degree(g))
function Padé(f, g, k::Integer)
    @assert degree(f) <= degree(g)
    r, t, s = fastconstrainedEEA(g, f, k)
    r, s, t
end

function slowgcd(g, f)
    R = parent(g)  # = K[x]
    U = (one(R), zero(R), f)  # = (1, 0, f)
    V = (zero(R), one(R), g)  # = (0, 1, g)
    # in Nemo, degree(0) is -1
    while !iszero(V[3])
        q = div(U[3], V[3])
        T = U .- q .* V
        U = V
        V = T
    end
    U[3]
end

# given (polynomials) g and f (|f| >= |g|),
# computes and returns a single row from the EEA algorithm (r, t, s), 
# such that r = t*g + s*f, |r| < k, where |r| is the maximal possible
function constrainedEEA(g, f, k::Integer)
    @assert degree(g) >= degree(f)
    R = parent(g)  # = K[x]
    U = (one(R), zero(R), f)  # = (1, 0, f)
    V = (zero(R), one(R), g)  # = (0, 1, g)
    # in Nemo, degree(0) is -1
    while degree(V[3]) > k
        q = div(U[3], V[3])
        T = U .- q .* V
        U = V
        V = T
    end
    (V[3], V[2], V[1])
end
