### ParamPanPam.jl - parametric Groebner bases.

### What

Groebner bases over $\mathbb{K}(a)[x]$, where $\mathbb{K}$ is either $\mathbb{Z}_p$ or $\mathbb{Q}$, and $a$ are treated as transcendental numbers.

Install `ParamPanPam.jl` with the following command in Julia

```julia
import Pkg; Pkg.add(url="https://github.com/sumiya11/ParamPanPam.jl")
```

### How

```julia
using Nemo, ParamPanPam

Ra, (a, b) = PolynomialRing(QQ, [:a, :b])
Rb, (x, y) = PolynomialRing(Rb, [:x, :y])

F = F = [b*x^2 + a//b, x*y - 5a];

ParamPanPam.paramgb(F)
[ Info: Given 2 polynomials in K(y)[x]
[ Info: Variables: AbstractAlgebra.Generic.MPoly{AbstractAlgebra.Generic.Frac{fmpq_mpoly}}[x, y]  
[ Info: Parameters: fmpq_mpoly[a, b]
┌ Info: Lifting to K[y][x]
│   fractionfreepolys =
│    2-element Vector{AbstractAlgebra.Generic.MPoly{fmpq_mpoly}}:
│     b^2*x^2 + a
└     x*y - 5*a
[ Info: Specializing at 1 + 2 random points to guess the shape..
┌ Info: The shape of the basis is: 2 polynomials with monomials
│   shapegb.shape =
│    2-element Vector{Vector{fmpq_mpoly}}:
│     [y^2, 1]
└     [x, y]
[ Info: Specializing at random points to guess the exponents in the coefficients..
┌ Info: The exponents in the coefficients
│   degrees =
│    2-element Vector{Vector{Tuple{Int64, Int64}}}:
│     [(0, 0), (3, 0)]
└     [(0, 0), (0, 2)]
[ Info: Initializing interpolation routines..
[ Info: Interpolation is bound by degrees 3, 3
[ Info: 1 points used..
[ Info: 2 points used..
[ Info: 4 points used..
[ Info: 8 points used..
┌ Warning: Coefficient basis[1][1] is interpolated!
│   P // Q = 1
└ @ Main.ParamPanPam c:\data\projects\ParamPanPam.jl\src\gb.jl:150
┌ Warning: Coefficient basis[2][1] is interpolated!
│   P // Q = 1
└ @ Main.ParamPanPam c:\data\projects\ParamPanPam.jl\src\gb.jl:150
[ Info: 16 points used..
┌ Warning: Coefficient basis[1][2] is interpolated!
│   P // Q = 25*a*b^2
└ @ Main.ParamPanPam c:\data\projects\ParamPanPam.jl\src\gb.jl:150
┌ Warning: Coefficient basis[2][2] is interpolated!
│   P // Q = 1//5//b^2
└ @ Main.ParamPanPam c:\data\projects\ParamPanPam.jl\src\gb.jl:150
[ Info: Success! 24 points used.
2-element Vector{AbstractAlgebra.Generic.MPoly{AbstractAlgebra.Generic.Frac{fmpq_mpoly}}}:
 y^2 + 25*a*b^2
 x + 1//5//b^2*y
```
