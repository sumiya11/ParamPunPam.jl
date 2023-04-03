#=
    Cuyt-Lee multivariate rational function interpolation.
=#

mutable struct FasterCuytLee{Ring, UnivRing}
    ring::Ring
    # the total degrees of the numerator/denominator
    Nd::Int
    Dd::Int
    # the number of terms in the numerator/denominator
    Nt::Int
    Dt::Int
    # dense interpolator for univariate functions
    cauchy::FasterCauchy{UnivRing}
    # polynomial interpolators for the numerator/denominator
    Ni::PrimesBenOrTiwari{Ring}
    Di::PrimesBenOrTiwari{Ring}

    T
    shift
    dilation
    ω
    ωs
    ξij_T
    points

    J

    function FasterCuytLee(
        ring::Ring,
        Nd::Int, Dd::Int,
        Nds::Vector{<:Integer}, Dds::Vector{<:Integer},
        Nt::Int, Dt::Int
        ) where {Ring}
        @assert Nd >= 0 && Dd >= 0
        @assert Nt >= 0 && Dt >= 0
        @assert all(>=(0), Nds) && all(>=(0), Dds)
        n = nvars(ring)
        @assert length(Dds) == length(Nds) == n
        K = base_ring(ring)
        Runiv, _ = Nemo.PolynomialRing(K, "u")
        cauchy = FasterCauchy(Runiv, Nd, Dd)
        Ni = PrimesBenOrTiwari(ring, Nt, Nd)
        Di = PrimesBenOrTiwari(ring, Dt, Dd)

        shift = random_point(ring)
        dilation = random_point(ring)
        ω = startingpoint(Ni)

        T = max(Nt, Dt)

        new{Ring,typeof(Runiv)}(
            ring, 
            Nd, Dd, 
            Nt, Dt, 
            cauchy, 
            Ni, Di,
            T, shift, dilation, ω, 
            map(i -> ω .^ i, 0:2T-1), 
            Vector{Vector{elem_type(K)}}(undef, 2T),
            Vector{Vector{elem_type(K)}}(undef, 2T*(Nd + Dd + 2)),
            -1
        )
    end
end

function get_evaluation_points!(cl::FasterCuytLee)
    if cl.J == -1
        cl.J = 0
    else
        cl.J = 2*cl.T
        cl.Nt *= 2
        cl.Dt *= 2
        cl.T = max(cl.Nt, cl.Dt)
        ringhom = cl.Ni.ring
        cl.Ni = PrimesBenOrTiwari(ringhom, cl.Nt, cl.Nd)
        cl.Di = PrimesBenOrTiwari(ringhom, cl.Dt, cl.Dd)
        resize!(cl.ξij_T, 2cl.T)
        resize!(cl.points, 2cl.T*(cl.Nd + cl.Dd + 2))
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
    cl.ωs = ωs = map(i -> ω .^ i, 0:2T-1)
    # ... and the sequence itself (the first degree is 0)
    J = cl.J

    used = Dict{elem_type(K), Bool}()
    ξij = distinct_points(K, Nd + Dd + 2)
    ωξij = Vector{Vector{elem_type(K)}}(undef, Nd + Dd + 2)
    ωξij0 = Vector{elem_type(K)}(undef, Nd + Dd + 2)

    @inbounds for i in J:2T-1
        ωi = ωs[i + 1]
        used[ξij[1]] = true
        while haskey(used, ξij[1])
            ξij[1] = random_point(K)
        end
        # The cycle below is (D + 2)*4n*log(q)
        for j in 1:Nd + Dd + 2
            !isassigned(ωξij, j) && (ωξij[j] = zeros(K, length(ω)))
            for nj in 1:length(ω)
                ωξij[j][nj] = ωi[nj]
            end
            ξ = ξij[j]
            for nj in 1:length(ω)
                ωξij[j][nj] = ωξij[j][nj]*ξ + shift[nj]
            end
        end
        cl.ξij_T[i + 1] = [x for x in ξij]
        cl.points[(i)*(Nd + Dd + 2) + 1 : (i + 1)*(Nd + Dd + 2)] .= deepcopy(ωξij)
    end
    cl.points
end

function interpolate!(cl::FasterCuytLee, evaluations::Vector{OO}) where {OO}
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
    P_coeffs = [Vector{elem_type(K)}(undef, 2T) for _ in 0:Nd]
    Q_coeffs  = [Vector{elem_type(K)}(undef, 2T) for _ in 0:Dd]
    P_interpolated = Vector{elem_type(R)}(undef, Nd+1)
    Q_interpolated = Vector{elem_type(R)}(undef, Dd+1)
    P_higher_degrees_contribution = [zeros(K, Nd+1) for _ in 0:2T-1]
    Q_higher_degrees_contribution = [zeros(K, Dd+1) for _ in 0:2T-1]
    for i in 0:2T - 1
        fij = evaluations[(i)*(Nd + Dd + 2) + 1 : (i + 1)*(Nd + Dd + 2)]
        # interpolate the numerator and the denominator densely.
        # M(D)logD
        P, Q = interpolate!(cauchy, ξij_T[i + 1], fij)
        @assert isone(trailing_coefficient(Q))
        for (poly, cfs) in ((P, P_coeffs), (Q, Q_coeffs))
            for jj in 1:length(cfs)
                cfs[jj][i + 1] = coeff(poly, length(cfs) - jj)
            end
        end
    end

    for (cfs, interpolated, degreebound, higher_contributions) in (
                (P_coeffs, P_interpolated, Nd, P_higher_degrees_contribution),
                (Q_coeffs, Q_interpolated, Dd, Q_higher_degrees_contribution)
            )

        # the index of the coefficient being intepolated            
        for idx in 1:degreebound + 1
            # the degree of the coefficient being intepolated
            deg = degreebound - idx + 1
            # evaluations at points ωs
            y_points = cfs[idx]
            # account for the contribution of higher degree terms
            for i in 1:2T
                y_points[i] = y_points[i] - higher_contributions[i][deg+1]
            end
            # consider either all evaluation points
            interpolated[idx] = interpolate!(Ni, ωs[1:2T], y_points)
            
            # update the contributions of higher degree term expansions 
            # to the lower degree coefficients
            expansion = collect(terms(evaluate(interpolated[idx], xs .+ shift)))
            # for each point..
            for ii in 1:2T
                # for each degree smaller than the current..
                for dd in 0:deg-1
                    td = filter(t -> total_degree(t) == dd, expansion)
                    tdii = sum(td, init=zero(R))
                    higher_contributions[ii][dd+1] += evaluate(tdii, ωs[ii])
                end
            end

        end
    end
    P = sum(P_interpolated)
    Q = sum(Q_interpolated)
    # undilated = xs .* map(inv, dilation) 
    # P = evaluate(P, undilated)
    # Q = evaluate(Q, undilated)
    normalization_factor = trailing_coefficient(Q)
    if !iszero(normalization_factor)
        P = map_coefficients(c -> div(c, normalization_factor), P)
        Q = map_coefficients(c -> div(c, normalization_factor), Q)
    end
    P, Q
end

function interpolate!(cl::FasterCuytLee, blackbox)
    R = cl.ring
    xs = Nemo.gens(R)
    K = base_ring(R)
    N, D = cl.N, cl.D
    uri = cl.univariate_rational_interpolator
    mpi = cl.multivariate_poly_interpolator
    # arrays that store arrays of evaluated coefficients
    # of numerator and denominator
    P_coeffs = [Vector{elem_type(K)}(undef, 0) for _ in 0:N]
    Q_coeffs  = [Vector{elem_type(K)}(undef, 0) for _ in 0:D]
    P_higher_degrees_contribution = Vector{Vector{elem_type(K)}}(undef, 0)
    Q_higher_degrees_contribution = Vector{Vector{elem_type(K)}}(undef, 0)
    # arrays that stores interpolated coefficients
    # of numerator and denominator
    P_interpolated = Vector{elem_type(R)}(undef, N+1)
    Q_interpolated = Vector{elem_type(R)}(undef, D+1)
    # arrays that store interpolators
    P_interpolators = [copy(mpi) for _ in 0:N]
    Q_interpolators = [copy(mpi) for _ in 0:D]
    # indices of the coefficient currently being interpolated
    P_initial, Q_initial = Ref{Bool}(true), Ref{Bool}(true)
    P_r, Q_r = Ref{Int}(1), Ref{Int}(1)
    # random variable shift 
    shift = random_point(R)
    # multivariate polynomial interpolation points
    ωs = Vector{Vector{elem_type(K)}}(undef, 0)
    i = 0
    all_interpolated = false
    while !all_interpolated
        # next point for sparse polynomial interpolation
        ωi = next_point!(mpi, increment=true)
        push!(ωs, ωi)
        i += 1
        # random points for dense rational interpolation,
        # f(ξ*x1,ξ*x2,..., ξ*xn) for ξ in ξij
        ξij = [random_point(K) for _ in 0:N + D + 1]
        @assert allunique(ξij) 
        # "substitute" ω 
        # f(ξ*ωi1,ξ*ωi2,..., ξ*ωin) for ξ in ξij
        ωξij = [ωi .* ξ for ξ in ξij]
        # shift in each of the variables,
        # f(ξ*ωi1 + s1,ξ*ωi2 + s2,..., ξ*ωin + sn)
        ωξsij = [ωξ .+ shift for ωξ in ωξij]
        # evaluate the blackbox
        fij = map(blackbox, ωξsij)
        # interpolate the numerator and the denominator densely
        P, Q = interpolate!(uri, ξij, fij)
        # @info "" P Q
        @assert isone(trailing_coefficient(Q))
        # store coefficients of dense interpolation of P and Q
        # in P_coeffs and Q_coeffs respectively
        for (poly, cfs) in ((P, P_coeffs), (Q, Q_coeffs))
            for jj in 1:length(cfs)
                push!(cfs[jj], coeff(poly, length(cfs) - jj))
            end
        end
        # for the currect index r,        
        # interpolate the N - r (or D - r) coefficient,
        # simultaneously in numerator and denominator
        for (mark, cfs, interpolated, interpolators, r, initial, degreebound, higher_contributions) in (
                (:P, P_coeffs, P_interpolated, P_interpolators, P_r, P_initial, N, P_higher_degrees_contribution),
                (:Q, Q_coeffs, Q_interpolated, Q_interpolators, Q_r, Q_initial, D, Q_higher_degrees_contribution)
            )
            push!(higher_contributions, zeros(K, degreebound + 1))
            success = true
            while success && r[] <= degreebound + 1
                # the index of the coefficient being intepolated
                idx = r[]
                # the degree of the coefficient being intepolated
                deg = degreebound - idx + 1
                # the corresponding intepolator object
                interpolator = interpolators[idx]
                # evaluations at points ωs
                y_points = cfs[idx]
                # account for the contribution of higher degree terms
                for i in 1:length(y_points)
                    y_points[i] = y_points[i] - higher_contributions[i][deg+1]
                end
                # consider either all or only the last evaluation point
                # (one point in case all previous points were
                # considered on the previous iteration of the outer while loop)
                !initial[] && (y_points = y_points[end:end])
                initial[] = false
                for y_point in y_points
                    _ = next_point!(interpolator)
                    success, f = next!(interpolator, y_point, y_point)
                    interpolated[idx] = f
                end
                # @warn "$mark" interpolated[idx]
                # if the coefficient is successfuly interpolated
                if success
                    # update the contributions of higher degree term expansions 
                    # to the lower degree coefficients
                    expansion = collect(terms(evaluate(interpolated[idx], xs .+ shift)))
                    # for each previously visited point..
                    for ii in 1:i
                        # for each degree smaller than the current..
                        for dd in 0:deg-1
                            td = filter(t -> total_degree(t) == dd, expansion)
                            tdii = sum(td, init=zero(R))
                            higher_contributions[ii][dd+1] += evaluate(tdii, ωs[ii])
                        end
                    end
                    r[] += 1
                    initial[] = true
                end
            end
        end
        all_interpolated = P_r[] + Q_r[] == N + D + 4
        if i > 2^7
            throw(ErrorException("Something bad happened in CuytLee"))
        end
    end
    # ans = P//Q
    P = sum(P_interpolated)
    Q = sum(Q_interpolated)
    normalization_factor = trailing_coefficient(Q)
    P = map_coefficients(c -> div(c, normalization_factor), P)
    Q = map_coefficients(c -> div(c, normalization_factor), Q)
    P, Q
end
