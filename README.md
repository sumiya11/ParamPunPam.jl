### ParamPanPam.jl - parametric Groebner bases.

### What

Groebner bases over $\mathbb{K}(a)[x]$, where $\mathbb{K}$ is either $\mathbb{Z}_p$ or $\mathbb{Q}$, and $a$ are treated as transcendental numbers.

You can install `ParamPanPam.jl` with the following command in Julia

```julia
import Pkg; Pkg.add(url="https://github.com/sumiya11/ParamPanPam.jl")
```

### How

```julia
using ParamPanPam
using Nemo

Ra, (a,b,c) = PolynomialRing(QQ, ["a","b","c"])
Rx, (x, y, z) = PolynomialRing(Nemo.FractionField(Ra), ["x", "y", "z"], ordering=:degrevlex)

F = [
    a*x^2 + b^2*x + (a + 1),
    x*y + b*y*z + 1//(a*b),
    x*z + c*z + b
];

ParamPanPam.paramgb(F)
[ Info: Given 3 polynomials in K(y)[x]
[ Info: Variables: AbstractAlgebra.Generic.MPoly{AbstractAlgebra.Generic.Frac{fmpq_mpoly}}[x, y, z]
[ Info: Parameters: fmpq_mpoly[a, b, c]
┌ Info: Lifting to K[y][x]
│   fractionfreepolys =
│    3-element Vector{AbstractAlgebra.Generic.MPoly{fmpq_mpoly}}:
│     a*x^2 + b^2*x + a + 1
│     a*b*x*y + a*b^2*y*z + 1
└     x*z + c*z + b
[ Info: Specializing at 1 + 2 random points to guess the shape..
┌ Info: The shape of the basis is: 3 polynomials with monomials
│   shapegb.shape =
│    3-element Vector{Vector{fmpq_mpoly}}:
│     [y, z, 1]
│     [x, z, 1]
└     [z^2, z, 1]
[ Info: Specializing at random points to guess the exponents in the coefficients..
[ Info: 1 points used..
[ Info: 2 points used..
[ Info: 4 points used..
[ Info: 8 points used..
[ Info: 16 points used..
[ Info: Success! 21 points used.
┌ Info: The exponents in the coefficients
│   degrees =
│    3-element Vector{Vector{Tuple{Int64, Int64}}}:
│     [(0, 0), (6, 9), (4, 7)]
│     [(0, 0), (3, 2), (2, 1)]
└     [(0, 0), (3, 3), (3, 3)]
[ Info: Initializing interpolation routines..
[ Info: Interpolation is bound by degrees 6, 9
[ Info: 1 points used..
[ Info: 2 points used..
[ Info: 4 points used..
[ Info: 8 points used..
[ Info: 16 points used..
[ Info: 32 points used..
┌ Info: Coefficient basis[1][1] is interpolated!
└   P // Q = 1
┌ Info: Coefficient basis[2][1] is interpolated!
└   P // Q = 1
┌ Info: Coefficient basis[3][1] is interpolated!
└   P // Q = 1
[ Info: 64 points used..
┌ Info: Coefficient basis[2][3] is interpolated!
└   P // Q = (-a*c + b^2)//a
[ Info: 128 points used..
┌ Info: Coefficient basis[2][2] is interpolated!
└   P // Q = (-a*c^2 - a + b^2*c - 1)//(a*b)
┌ Info: Coefficient basis[3][2] is interpolated!
└   P // Q = (2*a*b*c - b^3)//(a*c^2 + a - b^2*c + 1)
┌ Info: Coefficient basis[3][3] is interpolated!
└   P // Q = (a*b^2)//(a*c^2 + a - b^2*c + 1)
[ Info: 256 points used..
┌ Info: Coefficient basis[1][3] is interpolated!
└   P // Q = (-2*a*b^2*c - a*c^3 - a*c + b^4 + b^2*c^2 - c)//(a^2*b^5 + 2*a^2*b^3 + a^2*b*c^2 + a^2*b + a*b^5*c - a*b^3*c + 2*a*b^3 + a*b*c^2 + 2*a*b - b^7 - b^3*c + b)
┌ Info: Coefficient basis[1][2] is interpolated!
└   P // Q = (-a^2*b^2*c^2 - a^2*b^2 - a^2*c^4 - 2*a^2*c^2 - a^2 + a*b^4*c + 2*a*b^2*c^3 + 2*a*b^2*c - a*b^2 - 2*a*c^2 - 2*a - b^4*c^2 + 2*b^2*c - 1)//(a^3*b^6 + 2*a^3*b^4 + a^3*b^2*c^2 + a^3*b^2 + a^2*b^6*c - a^2*b^4*c + 2*a^2*b^4 
+ a^2*b^2*c^2 + 2*a^2*b^2 - a*b^8 - a*b^4*c + a*b^2)
[ Info: Success! 493 points used.
```

### Why




