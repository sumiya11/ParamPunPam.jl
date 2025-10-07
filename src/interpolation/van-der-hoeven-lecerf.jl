# van-der-Hoeven & Lecerf rational multivariate interpolation.

# Needs to know the total degrees, partial degrees, and the number of terms
# of numerator and denominator of the interpolant.
mutable struct VanDerHoevenLecerf{Ring, UnivRing, FiniteFieldElem}
    # multivariate polynomial ring
    ring::Ring
    # the total degrees of the numerator/denominator
    Nd::Int
    Dd::Int
    # vectors of partial degrees in the numerator/denominator
    Nds::Vector{Int}
    Dds::Vector{Int}
    # the number of terms in the numerator/denominator
    Nt::Int
    Dt::Int
    # dense interpolator for univariate functions
    cauchy::CauchyInterpolator{UnivRing}
    # polynomial interpolators for the numerator/denominator
    Ni::PrimesBenOrTiwari{Ring}
    Di::PrimesBenOrTiwari{Ring}

    T::Int
    shift::Vector{FiniteFieldElem}
    dilation::Vector{FiniteFieldElem}
    ω::Vector{FiniteFieldElem}
    ωs::Vector{Vector{FiniteFieldElem}}
    ξij_T::Vector{Vector{FiniteFieldElem}}
    ωξij0_T::Vector{Vector{FiniteFieldElem}}
    points::Vector{Vector{FiniteFieldElem}}

    J::Int

    function VanDerHoevenLecerf(
        ring::Ring,
        Nd::Int,
        Dd::Int,
        Nds::Vector{<:Integer},
        Dds::Vector{<:Integer},
        Nt::Int,
        Dt::Int
    ) where {Ring}
        @assert Nd >= 0 && Dd >= 0
        @assert Nt >= 0 && Dt >= 0
        @assert all(>=(0), Nds) && all(>=(0), Dds)
        n = nvars(ring)
        @assert length(Dds) == length(Nds) == n
        K = base_ring(ring)
        Runiv, _ = Nemo.polynomial_ring(K, "u")
        cauchy = CauchyInterpolator(Runiv, Nd, Dd)
        ringhom = homogenize(ring)
        Nds = vcat(Nd, Nds)
        Dds = vcat(Dd, Dds)
        # (!) we take the maximum for each of the partial degrees,
        # in order to be able use the same points for numerator/denominator
        NDds = map(maximum, zip(Nds, Dds))
        Ni = PrimesBenOrTiwari(ringhom, Nt, Nd)
        Di = PrimesBenOrTiwari(ringhom, Dt, Dd)
        shift = random_point(ringhom)
        dilation = random_point(ringhom)
        ω = startingpoint(Ni)

        T = max(Nt, Dt)

        new{Ring, typeof(Runiv), elem_type(K)}(
            ring,
            Nd,
            Dd,
            Nds,
            Dds,
            Nt,
            Dt,
            cauchy,
            Ni,
            Di,
            T,
            shift,
            dilation,
            ω,
            map(i -> ω .^ i, 0:(2T - 1)),
            Vector{Vector{elem_type(K)}}(undef, 2T),
            Vector{Vector{elem_type(K)}}(undef, 2T),
            Vector{Vector{elem_type(K)}}(undef, 2T * (Nd + Dd + 2)),
            -1
        )
    end
end

function get_evaluation_points!(vdhl::VanDerHoevenLecerf)
    if vdhl.J == -1
        vdhl.J = 0
    else
        vdhl.J = 2 * vdhl.T
        vdhl.Nt *= 2
        vdhl.Dt *= 2
        vdhl.T = max(vdhl.Nt, vdhl.Dt)
        ringhom = vdhl.Ni.ring
        vdhl.Ni = PrimesBenOrTiwari(ringhom, vdhl.Nt, vdhl.Nd)
        vdhl.Di = PrimesBenOrTiwari(ringhom, vdhl.Dt, vdhl.Dd)
        resize!(vdhl.ξij_T, 2vdhl.T)
        resize!(vdhl.ωξij0_T, 2vdhl.T)
        resize!(vdhl.points, 2vdhl.T * (vdhl.Nd + vdhl.Dd + 2))
    end

    R = vdhl.ring
    K = base_ring(R)
    Nd, Dd = vdhl.Nd, vdhl.Dd
    Nt, Dt = vdhl.Nt, vdhl.Dt
    Ni, Di = vdhl.Ni, vdhl.Di

    T = vdhl.T
    # Polynomial ring K[x0,x1,x2...xn]
    Rhom = Ni.ring
    # We will substitute points, such that
    # point[j]*dilation[j] + shift[j]
    shift = vdhl.shift
    dilation = vdhl.dilation
    # The starting point in the geometric sequence...
    ω = vdhl.ω
    vdhl.ωs = ωs = map(i -> ω .^ i, 0:(2T - 1))
    # ... and the sequence itself (the first degree is 0)
    J = vdhl.J
    totaldeg = Nd + Dd + 2

    used = Dict{elem_type(K), Bool}()
    ξij = distinct_nonzero_points(K, totaldeg)
    ωξij = Vector{Vector{elem_type(K)}}(undef, totaldeg)
    ωξij0 = Vector{elem_type(K)}(undef, totaldeg)

    # This cycle below is 
    # T*((D + 2)*4n*log(q) + D*L + 2*D*log(q) + M(D)log(D)),
    # which is T*M(D)log(D) + O(T*n*D*log(q)) + T*D*L
    @inbounds for i in J:(2T - 1)
        ωi = ωs[i + 1]
        ω0 = ωi[1]
        used[ξij[1]] = true
        while haskey(used, ξij[1])
            ξij[1] = random_point(K)
        end
        # The cycle below is (D + 2)*4n*log(q)
        for j in 1:(totaldeg)
            !isassigned(ωξij, j) && (ωξij[j] = zeros(K, length(ω) - 1))
            for nj in 2:length(ω)
                ωξij[j][nj - 1] = ωi[nj]
            end
            ξ = ξij[j]
            ωξij0[j] = ω0 * ξ * dilation[1] + shift[1]
            for nj in 2:length(ω)
                ωξij[j][nj - 1] = ωξij[j][nj - 1] * ξ * dilation[nj] + shift[nj]
                ωξij[j][nj - 1] = ωξij[j][nj - 1] // ωξij0[j]
            end
        end
        vdhl.ξij_T[i + 1] = [x for x in ξij]
        vdhl.ωξij0_T[i + 1] = [x for x in ωξij0]
        for j in ((i) * totaldeg + 1):((i + 1) * totaldeg)
            vdhl.points[j] = Vector{elem_type(K)}(undef, length(ω) - 1)
            for k in 1:(length(ω) - 1)
                vdhl.points[j][k] = ωξij[j - ((i) * totaldeg + 1) + 1][k]
            end
        end
    end

    vdhl.points
end

function interpolate!(vdhl::VanDerHoevenLecerf, evaluations::Vector{FiniteFieldElem}) where {FiniteFieldElem}
    T = vdhl.T
    R = vdhl.ring
    K = base_ring(R)
    Nd, Dd = vdhl.Nd, vdhl.Dd
    Nt, Dt = vdhl.Nt, vdhl.Dt
    Ni, Di = vdhl.Ni, vdhl.Di
    cauchy = vdhl.cauchy
    ωξij0_T = vdhl.ωξij0_T
    ξij_T = vdhl.ξij_T
    ωs = vdhl.ωs
    dilation = vdhl.dilation
    Rhom = Ni.ring
    Nys = Vector{elem_type(K)}(undef, 2T)
    Dys = Vector{elem_type(K)}(undef, 2T)
    @inbounds for i in 0:(2T - 1)
        fij = evaluations[((i) * (Nd + Dd + 2) + 1):((i + 1) * (Nd + Dd + 2))]
        fij = map(cξ -> cξ[1] * cξ[2]^(Nd - Dd), zip(fij, ωξij0_T[i + 1]))
        # interpolate the numerator and the denominator densely.
        # M(D)logD
        N, D = interpolate!(cauchy, ξij_T[i + 1], fij)
        @assert isone(trailing_coefficient(D))
        Nys[i + 1] = leading_coefficient(N)
        Dys[i + 1] = leading_coefficient(D)
    end
    # Interpolate the leading coefficients in the numerator and denominator
    # M(T)logT
    success_num, num = interpolate!(Ni, ωs[1:(2 * Nt)], Nys[1:(2 * Nt)])
    success_den, den = interpolate!(Di, ωs[1:(2 * Dt)], Dys[1:(2 * Dt)])
    success = success_num && success_den
    if !success
        return success, one(R), one(R)
    end
    # backward dilation,
    # substitute (x0,x1,x2,...xn) = (inv(d0)x0,inv(d1)x1,...,inv(dn)xn)
    # n*log(q), and T
    undilated = gens(Rhom) .* map(inv, dilation)
    num = evaluate(num, undilated)
    den = evaluate(den, undilated)
    # dehomogenization,
    # substitute (x0,x1,x2...,xn) = (1,x1,x2...,xn),
    xs0 = vcat([one(R)], gens(R))
    num = evaluate(num, xs0)
    den = evaluate(den, xs0)
    # normalize by the trailing_coefficient,
    # T*n*log(q)
    normalization_factor = iszero(den) ? zero(R) : trailing_coefficient(den)
    if !iszero(normalization_factor)
        num = map_coefficients(c -> div(c, normalization_factor), num, parent=R)
        den = map_coefficients(c -> div(c, normalization_factor), den, parent=R)
    end
    success, num, den
end
