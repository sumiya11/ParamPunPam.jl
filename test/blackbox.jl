using Nemo

@testset "Blackbox" begin
    Ra, (a,) = PolynomialRing(Nemo.QQ, ["a"], ordering=:degrevlex)
    Rx, (x, y, z) = PolynomialRing(Nemo.FractionField(Ra), ["x", "y", "z"], ordering=:degrevlex)

    s = [(1 // a)*x + (1 // (a + 1))*y, x + a//(a + 1)]
    bb = ParamPunPam.BasicBlackboxIdeal(s)
    @test parent(bb) == first(PolynomialRing(Ra, ["x", "y", "z"], ordering=:degrevlex))
    @test ParamPunPam.parent_params(bb) == Ra
    @test base_ring(bb) == Nemo.QQ
    K = Nemo.GF(2^31-1)
    ParamPunPam.reduce_mod_p!(bb, K)
    point = map(K, [1])
    # @test ParamPunPam.specialize_mod_p(bb, point) == [2x + y, 2x + 1]
end
