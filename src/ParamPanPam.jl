module ParamPanPam

using Nemo
using Groebner
import Primes

import Pkg
if !haskey(Pkg.installed(), "ExactSparseInterpolations")
    Pkg.add(url="https://github.com/sumiya11/ExactSparseInterpolations.jl")
end
import ExactSparseInterpolations

# Rational reconstructions
include("reconstruct.jl")

# Stuff
include("long-forgotten-truth.jl")

# The main algorithm
# include("gb.jl")

include("buchberger.jl")

export paramgb

end
