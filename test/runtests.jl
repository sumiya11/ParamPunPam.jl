using Test
using TestSetExtensions

import AbstractAlgebra
using Nemo

using ParamPunPam

@testset "All tests" verbose=true begin
    @includetests ["gb"]
    @includetests ["noon"]
end
