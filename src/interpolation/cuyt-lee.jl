#=
    Cuyt-Lee multivariate rational function interpolation.
=#

mutable struct CuytLee{Ring, I1<:AbstractRationalInterpolator, I2<:AbstractPolynomialInterpolator} <: AbstractInterpolator
    ring::Ring
    N::Int
    D::Int
    univariate_rational_interpolator::I1
    multivariate_poly_interpolator::I2
end

function CuytLee(ring, N::Integer, D::Integer; 
        univariate_rational_interpolator=DirectSolveRational(univariatize(AbstractAlgebra.PolyRing, ring), N, D),
        multivariate_poly_interpolator=BenOrTiwari(ring))
    CuytLee(ring, N, D, 
        univariate_rational_interpolator, 
        multivariate_poly_interpolator)
end

function interpolate!(cl::CuytLee, blackbox)
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
