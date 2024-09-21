using Nemo

@testset "Blackbox" begin
    Ra, (a,) = polynomial_ring(Nemo.QQ, ["a"], internal_ordering=:degrevlex)
    Rx, (x, y, z) = polynomial_ring(Nemo.fraction_field(Ra), ["x", "y", "z"], internal_ordering=:degrevlex)

    s = [(1 // a) * x + (1 // (a + 1)) * y, x + a // (a + 1)]
    bb = ParamPunPam.BasicBlackboxIdeal(s)
    @test parent(bb) == first(polynomial_ring(Ra, ["x", "y", "z"], internal_ordering=:degrevlex))
    @test ParamPunPam.parent_params(bb) == Ra
    @test base_ring(bb) == Nemo.QQ
    K = Nemo.Native.GF(2^31 - 1)
    ParamPunPam.reduce_mod_p!(bb, K)
    point = map(K, [1])
    # @test ParamPunPam.specialize_mod_p(bb, point) == [2x + y, 2x + 1]
end
