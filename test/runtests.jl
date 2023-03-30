using Test
using TestSetExtensions

using Nemo

include("../src/ParamPunPam.jl")

@testset "All tests" verbose=true begin
    @includetests ["gb"]
end
