### ParamPunPam.jl - parametric Groebner bases.

### What is it?

ParamPunPam.jl computes Groebner bases in $\mathbb{Q}(\mathbf{a})[\mathbf{x}]$, where $\mathbf{a}$ are treated as transcendental numbers.

ParamPunPam.jl is primarily designed to work in cases where the parametric coefficients of the output basis are sparse.

You can install ParamPunPam.jl with the following command in Julia:

```julia
import Pkg; Pkg.add("ParamPunPam.jl")
```

### How to use it?

ParamPunPam.jl works in combination with polynomials from Nemo.jl and provides a single command `paramgb`.
See the following example:

```julia
using ParamPunPam, Nemo

R_param, (a, b) = QQ["a", "b"]
R, (x, y, z) = FractionField(R_param)["x", "y", "z"]

F = [
    x^2 + x + (a + 1),
    x*y + b*y*z + 1//(a*b),
    x*z + z + b
]

paramgb(F)

# returns
3-element Vector{AbstractAlgebra.Generic.MPoly{AbstractAlgebra.Generic.Frac{fmpq_mpoly}}}:
 z^2 + b//(a + 1)*z + b^2//(a + 1)
 y + (-a - 1)//(a^2*b^2 + a*b^4 + a*b^2)*z - 1//(a^2*b + a*b^3 + a*b)
 x + (-a - 1)//b*z
```
