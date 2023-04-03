using Test
using TestSetExtensions

import AbstractAlgebra
using Nemo

include("../src/ParamPunPam.jl")

@testset "All tests" verbose=true begin
    @includetests ["gb"]
    @includetests ["noon"]
end
