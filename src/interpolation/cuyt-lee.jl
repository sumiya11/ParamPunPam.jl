# Cuyt-Lee multivariate rational function interpolation.

mutable struct CuytLee{Ring, UnivRing, FiniteFieldElem}
    ring::Ring
    # the total degrees of the numerator/denominator
    Nd::Int
    Dd::Int
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
    points::Vector{Vector{FiniteFieldElem}}

    J::Int

    function CuytLee(
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
        Ni = PrimesBenOrTiwari(ring, Nt, Nd)
        Di = PrimesBenOrTiwari(ring, Dt, Dd)

        shift = random_point(ring)
        dilation = random_point(ring)
        # TODO!
        dilation = [one(K) for _ in 1:n]

        ω = startingpoint(Ni)

        T = max(Nt, Dt)

        new{Ring, typeof(Runiv), elem_type(K)}(
            ring,
            Nd,
            Dd,
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
            Vector{Vector{elem_type(K)}}(undef, 2T * (Nd + Dd + 2)),
            -1
        )
    end
end

function get_evaluation_points!(cl::CuytLee)
    if cl.J == -1
        cl.J = 0
    else
        cl.J = 2 * cl.T
        cl.Nt *= 2
        cl.Dt *= 2
        cl.T = max(cl.Nt, cl.Dt)
        ringhom = cl.Ni.ring
        cl.Ni = PrimesBenOrTiwari(ringhom, cl.Nt, cl.Nd)
        cl.Di = PrimesBenOrTiwari(ringhom, cl.Dt, cl.Dd)
        resize!(cl.ξij_T, 2cl.T)
        resize!(cl.points, 2cl.T * (cl.Nd + cl.Dd + 2))
    end

    R = cl.ring
    K = base_ring(R)
    Nd, Dd = cl.Nd, cl.Dd
    Nt, Dt = cl.Nt, cl.Dt
    Ni, Di = cl.Ni, cl.Di

    T = cl.T
    # Polynomial ring K[x0,x1,x2...xn]
    Rhom = Ni.ring
    # We will substitute points, such that
    # point[j]*dilation[j] + shift[j]
    shift = cl.shift
    dilation = cl.dilation
    # The starting point in the geometric sequence...
    ω = cl.ω
    cl.ωs = ωs = map(i -> ω .^ i, 0:(2T - 1))
    # ... and the sequence itself (the first degree is 0)
    J = cl.J
    totaldeg = Nd + Dd + 2

    @debug "Generating $(2T) evaluation points, $J..$(2T)"
    @debug """
    Params:
    shift = $shift
    dilation = $dilation
    sum of degrees = $totaldeg
    """

    used = Dict{elem_type(K), Bool}()
    ξij = distinct_nonzero_points(K, totaldeg)
    ωξij = Vector{Vector{elem_type(K)}}(undef, totaldeg)

    @inbounds for i in J:(2T - 1)
        ωi = ωs[i + 1]
        used[ξij[1]] = true
        while haskey(used, ξij[1])
            ξij[1] = random_point(K)
        end
        for j in 1:totaldeg
            !isassigned(ωξij, j) && (ωξij[j] = zeros(K, length(ω)))
            for nj in 1:length(ω)
                ωξij[j][nj] = ωi[nj]
            end
            ξ = ξij[j]
            for nj in 1:length(ω)
                ωξij[j][nj] = ωξij[j][nj] * ξ * dilation[nj] + shift[nj]
            end
        end
        cl.ξij_T[i + 1] = [x for x in ξij]
        for j in ((i) * totaldeg + 1):((i + 1) * totaldeg)
            cl.points[j] = Vector{elem_type(K)}(undef, length(ω))
            for k in 1:length(ω)
                cl.points[j][k] = ωξij[j - ((i) * totaldeg + 1) + 1][k]
            end
        end
    end

    @debug "Homogenizing variable" cl.ξij_T

    cl.points
end

function interpolate!(cl::CuytLee, evaluations::Vector{FiniteFieldElem}) where {FiniteFieldElem}
    T = cl.T
    R = cl.ring
    xs = gens(R)
    K = base_ring(R)
    Nd, Dd = cl.Nd, cl.Dd
    Nt, Dt = cl.Nt, cl.Dt
    Ni, Di = cl.Ni, cl.Di
    cauchy = cl.cauchy
    ξij_T = cl.ξij_T
    ωs = cl.ωs
    shift = cl.shift
    dilation = cl.dilation
    Rhom = Ni.ring

    totaldeg = Nd + Dd + 2
    P_coeffs = [Vector{elem_type(K)}(undef, 2T) for _ in 0:Nd]
    Q_coeffs = [Vector{elem_type(K)}(undef, 2T) for _ in 0:Dd]

    P_interpolated = Vector{elem_type(R)}(undef, Nd + 1)
    Q_interpolated = Vector{elem_type(R)}(undef, Dd + 1)
    P_higher_degrees_contribution = [zeros(K, Nd + 1) for _ in 0:(2T - 1)]
    Q_higher_degrees_contribution = [zeros(K, Dd + 1) for _ in 0:(2T - 1)]

    @inbounds for i in 0:(2T - 1)
        fij = evaluations[((i) * totaldeg + 1):((i + 1) * totaldeg)]
        # interpolate the numerator and the denominator densely.
        # M(D)logD
        P, Q = interpolate!(cauchy, ξij_T[i + 1], fij)

        @debug "" i P Q

        @assert isone(trailing_coefficient(Q))
        for (poly, cfs) in ((P, P_coeffs), (Q, Q_coeffs))
            for jj in 1:length(cfs)
                cfs[jj][i + 1] = coeff(poly, length(cfs) - jj)
            end
        end
    end

    @debug "" P_coeffs Q_coeffs
    success = true

    for (cfs, interpolated, degreebound, higher_contributions, polynomial_interpolator, T_bound) in (
        (P_coeffs, P_interpolated, Nd, P_higher_degrees_contribution, Ni, Nt),
        (Q_coeffs, Q_interpolated, Dd, Q_higher_degrees_contribution, Di, Dt)
    )
        # the index of the coefficient being intepolated            
        @inbounds for idx in 1:(degreebound + 1)
            # the degree of the coefficient being intepolated
            deg = degreebound - idx + 1
            @debug """
            Interpolating coefficient at index $idx (out of $(degreebound + 1))
            Degree of that coefficient is $deg.
            """

            # evaluations at points ωs
            y_points = cfs[idx]

            @debug """
            Evaluation points: $(ωs[1:(2T)])
            Evaluation results: $(y_points)
            """

            @debug """
            Higher-order contributions:
            $higher_contributions
            """

            # account for the contribution of higher degree terms
            for i in 1:(2T)
                y_points[i] = y_points[i] - higher_contributions[i][deg + 1]
            end

            success_i, poly = interpolate!(polynomial_interpolator, ωs[1:(2 * T_bound)], y_points[1:(2 * T_bound)])
            interpolated[idx] = poly
            success = success_i && success

            @debug "Interpolated $(interpolated[idx])"
            # update the contributions of higher degree term expansions 
            # to the lower degree coefficients
            expansion = collect(terms(evaluate(interpolated[idx], (xs .* inv.(dilation)) .+ shift)))
            @debug """
            Expansion of interpolated coefficient:
            $expansion
            """

            # for each point..
            for i in 1:(2T)
                # for each degree smaller than the current..
                for dd in 0:(deg - 1)
                    td = filter(t -> total_degree(t) == dd, expansion)
                    tdi = sum(td, init=zero(R))
                    higher_contributions[i][dd + 1] += evaluate(tdi, ωs[i])
                    # higher_contributions[i][dd + 1] += cfs[idx][i]
                end
            end
        end
    end
    P = sum(P_interpolated)
    Q = sum(Q_interpolated)
    undilated = gens(Rhom) .* map(inv, dilation)
    P = evaluate(P, undilated)
    Q = evaluate(Q, undilated)
    normalization_factor = trailing_coefficient(Q)
    if !iszero(normalization_factor)
        P = map_coefficients(c -> div(c, normalization_factor), P, parent=R)
        Q = map_coefficients(c -> div(c, normalization_factor), Q, parent=R)
    end
    success, P, Q
end
