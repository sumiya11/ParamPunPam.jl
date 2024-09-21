
@testset "Regression: cancellation of leading terms" begin
    Ra, (a, b, c) = polynomial_ring(Nemo.QQ, ["a", "b", "c"])
    Rx, (x, y, z) = polynomial_ring(Nemo.fraction_field(Ra), ["x", "y", "z"], internal_ordering=:degrevlex)
    F = [x * y - y * z - a * b + a * c]
    @test ParamPunPam.paramgb(F) == F
    F = [x * y - y * z + (a * b)^3 + (a * c)^3 - (a * b)^2 + (a * c)^2 + a * b + a * c]
    @test ParamPunPam.paramgb(F) == F
end
