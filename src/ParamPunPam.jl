module ParamPunPam

import AbstractAlgebra

using Groebner
using Logging

using Nemo
import Primes

# For the fun!
using ProgressMeter

include("long-forgotten-truth.jl")
include("utils.jl")

include("interpolation/generic.jl")
include("interpolation/discrete-log.jl")
include("interpolation/fastgcd.jl")
include("interpolation/cauchy.jl")
include("interpolation/div-and-conq.jl")
include("interpolation/kron-ben-or-tiwari.jl")
include("interpolation/primes-ben-or-tiwari.jl")
include("interpolation/van-der-hoeven-lecerf.jl")
include("interpolation/cuyt-lee.jl")

include("groebner/blackbox.jl")
include("groebner/groebnerstate.jl")
include("groebner/reconstruct.jl")
include("groebner/modulartracker.jl")
include("groebner/paramgb.jl")

export paramgb
export AbstractBlackboxIdeal, BasicBlackboxIdeal

end
