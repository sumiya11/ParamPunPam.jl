module ParamPanPam

using Nemo
using Groebner

import Pkg
if !haskey(Pkg.installed(), "ExactSparseInterpolations")
    Pkg.add(url="https://github.com/sumiya11/ExactSparseInterpolations.jl")
end
using ExactSparseInterpolations

include("gb.jl")

export paramgb

end
