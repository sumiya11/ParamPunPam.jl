
@testset "GB over Q(a)" begin
    Ra, (a,) = PolynomialRing(QQ, ["a"])
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
    for (case, answer) in zip(cases, answers)
        G = ParamPanPam.paramgb(case)
        @warn "" case G
        @test G == answer 
    end

    Ra, (a,b,c) = PolynomialRing(QQ, ["a","b","c"])
    Rx, (x, y, z) = PolynomialRing(Nemo.FractionField(Ra), ["x", "y", "z"], ordering=:degrevlex)
    cases = [
        [x + (a + b + c)^3//(a*c*b)^2],
        [a*x^2 + b^2*x + (a + 1), 
         x*y + b*y*z + 1//(a*b), 
         x*z + c*z + b]
    ]
    answers = [
        [x + (a + b + c)^3//(a*c*b)^2],
        [y + (-a^2*b^2*c^2 - a^2*b^2 - 
        a^2*c^4 - 2*a^2*c^2 - a^2 + a*b^4*c + 2*a*b^2*c^3 + 2*a*b^2*c - a*b^2 - 2*a*c^2 - 2*a - b^4*c^2 + 2*b^2*c - 1)//(a^3*b^6 + 2*a^3*b^4 + a^3*b^2*c^2 + a^3*b^2 + a^2*b^6*c - a^2*b^4*c + 2*a^2*b^4 + a^2*b^2*c^2 + 2*a^2*b^2 - a*b^8 - a*b^4*c + a*b^2)*z + (-2*a*b^2*c - a*c^3 - a*c + b^4 + b^2*c^2 - c)//(a^2*b^5 + 2*a^2*b^3 + a^2*b*c^2 + a^2*b + a*b^5*c - a*b^3*c + 2*a*b^3 + a*b*c^2 + 2*a*b - b^7 - b^3*c + b), 
        x + (-a*c^2 - a + b^2*c - 1)//(a*b)*z + (-a*c + b^2)//a, 
        z^2 + (2*a*b*c - b^3)//(a*c^2 + a - b^2*c + 1)*z + (a*b^2)//(a*c^2 + a - b^2*c + 1)]
    ]
    for (case, answer) in zip(cases, answers)
        G = ParamPanPam.paramgb(case)
        @warn "" case G
        @test G == answer 
    end

end
