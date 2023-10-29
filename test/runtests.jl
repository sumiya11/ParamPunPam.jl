using Test
using TestSetExtensions

import AbstractAlgebra
using Nemo

using ParamPunPam

@testset "All tests" verbose = true begin
    @includetests ["discrete-log", "div-and-conq", "fastgcd"]
    @includetests ["ben-or-tiwari", "interpolators"]
    @includetests ["blackbox"]
    @includetests ["utils"]
    @includetests ["paramgb", "logging"]
    @includetests ["regressions"]
end
