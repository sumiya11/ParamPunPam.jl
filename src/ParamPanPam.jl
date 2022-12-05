module ParamPanPam

using Nemo
using Groebner

import Pkg; Pkg.add(url="https://github.com/sumiya11/ExactSparseInterpolations.jl")
import ExactSparseInterpolations

include("gb.jl")

export paramgb

end
