using Test
using TestSetExtensions

import AbstractAlgebra
using Nemo

using ParamPunPam

@testset "All tests" verbose=true begin
    @includetests ["fields", "discrete-log", "div-and-conq", "fastgcd"]
    @includetests ["ben-or-tiwari", "interpolators"]
    @includetests ["blackbox"]
    @includetests ["paramgb"]
    @includetests ["regressions"]
end
