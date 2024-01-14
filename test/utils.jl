@testset "Rational reconstruction" begin
    Rff, (a, b, c) = Nemo.Native.GF(2^31 - 1)["x", "y", "z"]
    Rqq, (aqq, bqq, cqq) = QQ["x", "y", "z"]
    cases = [
        (poly=Rff(1), ans=Rqq(1)),
        (poly=Rff(-100), ans=Rqq(-100)),
        (poly=2a + 3b - 5c, ans=2aqq + 3bqq - 5cqq),
        (poly=-10a + -10b - 10c + 99a^2, ans=-10aqq + -10bqq - 10cqq + 99aqq^2)
    ]
    for case in cases
        polyff = case.poly
        success, polyqq = ParamPunPam.rational_reconstruct_polynomial(Rqq, polyff)
        @test success && polyqq == case.ans
    end
end
