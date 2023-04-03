
function test_paramgb(cases, answers; kwargs...)
    for (case, answer) in zip(cases, answers)
        G = nothing
        # try
        #     G = ParamPunPam.paramgb(case)
        # catch e
        #     @warn "Exception: $e"
        #     rethrow(e)
        # end
        G = ParamPunPam.paramgb(case; kwargs...)
        @warn "" case G kwargs
        if haskey(kwargs, :up_to_degree)
            continue
        end 
        @test G == answer 
    end
end

@testset "GB over Q(a...)" begin
    Ra, (a,) = PolynomialRing(Nemo.QQ, ["a"], ordering=:degrevlex)
    Rx, (x, y, z) = PolynomialRing(Nemo.FractionField(Ra), ["x", "y", "z"], ordering=:degrevlex)
    cases = [
        [x, x + a],
        [(x + a)^2],
        [x + a^2, x*y + a^10],
        [x^2 - x + a*y^2 + a*z^2,
         a*x*y + a*y*z - y,
         x + a*y + a*z - 1]
    ]
    answers = [
        [Rx(1)],
        [(x + a)^2],
        [y - a^8, x + a^2],
        [x + a*y + a*z - 1,
         y*z + (a^2 + a)//(a^2 + 1)*z^2 - 1//(a^3 + a)*y - a//(a^2 + 1)*z,
         y^2 + (-a^2 + 1)//(a^2 + 1)*z^2 + (-a + 1)//(a^2 + 1)*y + (a - 1)//(a^2 + 1)*z,
         z^3 + (-3//2*a^5 + a^4 - a^3 + 1//2*a^2 - 1//2*a - 1//2)//(a^6 + 1//2*a^5 + a^4 + a^3 + 1//2*a)*z^2 + (1//2*a^3 - 1//2)//(a^6 + 1//2*a^5 + a^4 + a^3 + 1//2*a)*y + (1//2*a^4 - a^3 + 1//2*a^2 - 1//2*a + 1//2)//(a^6 + 1//2*a^5 + a^4 + a^3 + 1//2*a)*z]
    ]
    test_paramgb(cases, answers)
    test_paramgb(
        cases, answers,
        rational_interpolator=ParamPunPam.CuytLee()
    )
    test_paramgb(
        cases, answers,
        up_to_degree=(1, 1),
        rational_interpolator=ParamPunPam.CuytLee()
    )
    test_paramgb(
        cases, answers,
        up_to_degree=(4, 4),
    )

    Ra, (a1,a2,a3,a4,a5) = PolynomialRing(Nemo.QQ, ["a1","a2","a3","a4","a5"], ordering=:degrevlex)
    Rx, (x, y, z) = PolynomialRing(Nemo.FractionField(Ra), ["x", "y", "z"], ordering=:degrevlex)
    cases = [
        [(x - a1)*(y - a2)*(z - a3)*(x - a4)*(x - a5)],
    ]
    answers = [
        [(x - a1)*(y - a2)*(z - a3)*(x - a4)*(x - a5)],
    ]
    test_paramgb(cases, answers)

    Ra, (a,b,c) = PolynomialRing(QQ, ["a","b","c"])
    Rx, (x, y, z) = PolynomialRing(Nemo.FractionField(Ra), ["x", "y", "z"], ordering=:deglex)
    cases = [
        [x + (a + b + c)^3//(a*c*b)^2]
    ]
    answers = [
        [x + (a + b + c)^3//(a*c*b)^2]
    ]
    test_paramgb(cases, answers)
end
