using Nemo

@testset "Noon" begin
    Rparam, (a,b,c) = PolynomialRing(QQ, ["a","b","c"], ordering=:degrevlex)
    R, (x1, x2, x3) = PolynomialRing(FractionField(Rparam), ["x1","x2","x3"], ordering=:degrevlex)

    f = [
       a*x1*x2^2 + (a + c)*x1*x3^2 - b*x1 + a,
       a*x1^2*x2 + (a + b)*x2*x3^2 - b*x2 + a,
       a*x1^2*x3 + a*x2^2*x3 - b*x3 + a,
    ]
    
    ParamPunPam.paramgb(f)
end


