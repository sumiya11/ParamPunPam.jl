module ParamPanPam

using Nemo
using Groebner
import Primes

import Pkg
if !haskey(Pkg.installed(), "ExactSparseInterpolations")
    Pkg.add(url="https://github.com/sumiya11/ExactSparseInterpolations.jl")
end
import ExactSparseInterpolations

include("reconstruct.jl")
include("gb.jl")

export paramgb

end
