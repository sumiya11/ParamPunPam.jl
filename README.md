### ParamPunPam.jl - parametric Groebner bases.

[![Runtests](https://github.com/sumiya11/ParamPunPam.jl/actions/workflows/Runtests.yml/badge.svg)](https://github.com/sumiya11/ParamPunPam.jl/actions/workflows/Runtests.yml)

### What is it?

ParamPunPam.jl computes Groebner bases in $\mathbb{Q}(\mathbf{a})[\mathbf{x}]$, where $\mathbf{a}$ are treated as transcendental numbers.

ParamPunPam.jl is primarily designed to work in cases when the coefficients of the basis are sparse and of low degree.

You can install ParamPunPam.jl with the following command in Julia:

```julia
import Pkg; Pkg.add("ParamPunPam")
```

### How to use it?

ParamPunPam.jl provides a single function: `paramgb`.
This function works in combination with polynomials from Nemo.jl.
See the following example:

```julia
using ParamPunPam, Nemo

# Create polynomial rings
R_param, (a, b) = QQ["a", "b"]
R, (x, y, z) = fraction_field(R_param)["x", "y", "z"]

# Polynomial system
F = [
    x^2 + x + (a + 1),
    x*y + b*y*z + 1//(a*b),
    x*z + z + b
]

paramgb(F)

# Returns
3-element Vector{AbstractAlgebra.Generic.MPoly{AbstractAlgebra.Generic.FracFieldElem{fmpq_mpoly}}}:
 z^2 + b//(a + 1)*z + b^2//(a + 1)
 y + (-a - 1)//(a^2*b^2 + a*b^4 + a*b^2)*z - 1//(a^2*b + a*b^3 + a*b)
 x + (-a - 1)//b*z
```
