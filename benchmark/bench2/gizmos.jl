using StructuralIdentifiability
using RationalFunctionFields
using ParamPunPam
using Nemo
using Groebner

power_sum(x, k) = sum(x .^ k)
power_sums(x, k) = [power_sum(x, i) for i in 1:k]

function max_terms_of_coeffs(gb)
    isnothing(gb) && return nothing
    N, D = 0, 0
    for f in gb
        for c in collect(coefficients(f))
            N = max(N, length(numerator(c)))
            D = max(D, length(denominator(c)))
        end
    end
    N, D
end

function max_degree_of_coeffs(gb)
    isnothing(gb) && return nothing
    N, D = 0, 0
    for f in gb
        for c in collect(coefficients(f))
            N = max(N, total_degree(numerator(c)))
            D = max(D, total_degree(denominator(c)))
        end
    end
    N, D
end
