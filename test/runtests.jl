using Test
using TestSetExtensions

using Nemo

include("../src/ParamPanPam.jl")

@testset "All tests" verbose=true begin
    @includetests ["gb"]
end