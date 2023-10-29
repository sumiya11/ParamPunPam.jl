
@testset "Regression: cancellation of leading terms" begin
    Ra, (a, b, c) = PolynomialRing(Nemo.QQ, ["a", "b", "c"])
    Rx, (x, y, z) =
        PolynomialRing(Nemo.FractionField(Ra), ["x", "y", "z"], ordering=:degrevlex)
    F = [x * y - y * z - a * b + a * c]
    @test ParamPunPam.paramgb(F) == F
    F = [x * y - y * z + (a * b)^3 + (a * c)^3 - (a * b)^2 + (a * c)^2 + a * b + a * c]
    @test ParamPunPam.paramgb(F) == F
end
