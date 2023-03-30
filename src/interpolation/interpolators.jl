
mutable struct Interpolators{T, F}#, Iu, Im}
    # this are needed
    polys::Vector{T}
    # 
    univariate_rational_interpolator
    #
    multivariate_rational_interpolator
    # 
    statuses::Vector{Vector{Bool}}
    #
    degrees::Vector{Vector{Tuple{Int, Int}}}
    # 
    univ_η::Int
    # 
    mult_η::Int
    #
    univ_shift::F
    univ_dilate::F

    function Interpolators(polys::Vector{T}, modular; 
            univ_η::Integer=0, mult_η::Integer=0,
            univ_shift::Bool=true, univ_dilate::Bool=true
        ) where {T}
        @assert univ_η >= 0 && mult_η >= 0
        dilate = randluckypoint(modular)
        shift = randluckypoint(modular)
        new{T, typeof(dilate)}(
            polys, nothing, nothing, 
            Vector{Vector{Bool}}(), Vector{Vector{Tuple{Int, Int}}}(),
            univ_η, mult_η,
            dilate, shift
        )
    end
end

# function initialize_univariate_interpolator!(iss::Interpolators, univ_ring, shape)
#     iss.univariate_rational_interpolators = [
#         [
#             ExactSparseInterpolations.AdaptiveCauchy(univ_ring)
#             for _ in 1:length(shape[i])
#         ]
#         for i in 1:length(shape)
#     ]
#     iss.statuses = [
#         [false for _ in 1:length(shape[i])]
#         for i in 1:length(shape)
#     ]
#     iss.degrees = [
#         [(-1, -1) for _ in 1:length(shape[i])]
#         for i in 1:length(shape)
#     ]
#     npoints = 2
#     UnivariateDriver(iss, npoints, false)
# end

# function next_points!(ud::UnivariateDriver, state)
#     iss = ud.iss
#     one_of_a_kind = first(first(iss.univariate_rational_interpolators))
#     x_points = [
#         ExactSparseInterpolations.next_point!(one_of_a_kind)
#         for _ in 1:ud.npoints
#     ]
#     ud.npoints *= 2
#     map(x -> x * iss.univ_dilate + iss.univ_shift, x_points)
# end

# function next_evals!(ud::UnivariateDriver, xs, ys)
#     iss = ud.iss
#     for k in 1:length(ys)
#         basis = ys[k]
#         x_point = xs[k]
#         for i in 1:length(basis)
#             for j in 1:length(basis[i])
#                 status = iss.statuses[i][j]
#                 # if already interpolated
#                 status && continue
#                 interpolator = iss.univariate_rational_interpolators[i][j]
#                 y_point = basis[i][j]
#                 success, (P, Q) = ExactSparseInterpolations.next!(interpolator, x_point, y_point)
#                 iss.degrees[i][j] = (degree(P), degree(Q))
#                 iss.statuses[i][j] = success
#             end
#         end
#     end
#     ud.all_interpolated = all(all, iss.statuses)
#     nothing
# end
