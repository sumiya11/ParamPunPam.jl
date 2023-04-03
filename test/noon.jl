using Nemo

@testset "Noon" begin
    Rparam, (a,b,c) = PolynomialRing(QQ, ["a","b","c"])
    R, (x1, x2, x3) = PolynomialRing(FractionField(Rparam), ["x1","x2","x3"], ordering=:degrevlex)

    f = [
       a*x1*x2^2 + (a + c)*x1*x3^2 - b*x1 + a,
       a*x1^2*x2 + (a + b)*x2*x3^2 - b*x2 + a,
       a*x1^2*x3 + a*x2^2*x3 - b*x3 + a,
    ]
    
    ParamPunPam.paramgb(f)

    Rparam, (a1,a2,a3,a4,a5,a6) = PolynomialRing(QQ, ["a$i" for i in 1:6])
    R, (x1, x2, x3) = PolynomialRing(FractionField(Rparam), ["x1","x2","x3"], ordering=:degrevlex)

    f = [
       x1*x2^2 + (a1 + a2 + a3)*x1*x3^2 - (a2 + a5 + a6),
       x1^2*x2 + (a1 + a2)*x2*x3^2 - a2*x2 + a1,
       x1^2*x3 + a1*x2^2*x3 - a2*x3 + (a1 + a3 + a6),
    ]
    
    ParamPunPam.paramgb(f)

    Rparam, Ai = PolynomialRing(QQ, ["A$i" for i in 1:10])
    R, (x,y,z) = PolynomialRing(FractionField(Rparam), ["x","y","z"], ordering=:degrevlex)
    f = [
        sum(Ai)*x*y + (Ai[1] + Ai[10])//(sum(Ai)),
        y*z - sum(Ai),
        x^2 + Ai[5] // (sum(Ai))
    ]
    
    ParamPunPam.paramgb(f)
end


