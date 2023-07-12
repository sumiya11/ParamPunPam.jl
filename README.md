### ParamPunPam.jl - parametric Groebner bases.

### What is it?

Groebner bases in $\mathbb{Q}(\mathbf{a})[\mathbf{x}]$, where $\mathbf{a}$ are treated as transcendental numbers.

`ParamPunPam.jl` is primarily designed to work in cases where the coefficients of the output basis are sparse.

You can install `ParamPunPam.jl` with the following command in Julia:

```julia
import Pkg; Pkg.add(url="https://github.com/sumiya11/ParamPunPam.jl")
```

### How to use it?

See the following example:

```julia
using ParamPunPam, Nemo

Rparam, (a, b) = PolynomialRing(QQ, ["a", "b"])
R, (x, y, z) = PolynomialRing(FractionField(Rparam), ["x", "y", "z"], ordering=:degrevlex)

F = [
    x^2 + x + (a + 1),
    x*y + b*y*z + 1//(a*b),
    x*z + z + b
]

ParamPunPam.paramgb(BasicBlackboxIdeal(F))
[ Info: Constructing a BasicBlackboxIdeal from 3 input polynomials
┌ Info: Computing parametric GB
│   up_to_degree = (4611686018427387903, 4611686018427387903)
│   rational_interpolator = ParamPunPam.VanDerHoevenLecerf() 
└   polynomial_interpolator = :primes_benortiwari
[ Info: Given 3 functions in Rational Field(a, b)[x, y, z]
[ Info: Specializing at 1 + 2 random points to guess the basis shape..
[ Info: Reducing modulo Galois field with characteristic 4611686018427388039..
┌ Info: The shape of the basis is: 3 polynomials with monomials
│   state.shape =
│    3-element Vector{Vector{gfp_mpoly}}:
│     [y, z, 1]
│     [x, z]
└     [z^2, z, 1]
[ Info: Specializing at random points to guess the total degrees in parameters..
[ Info: Using 6 points..
[ Info: Using 10 points..
[ Info: Using 18 points..
[ Info: Using 34 points..
[ Info: Success! 34 points used.
┌ Info: The total degrees in the coefficients
│   state.param_degrees =
│    3-element Vector{Vector{Tuple{Int64, Int64}}}:
│     [(0, 0), (1, 5), (0, 4)]
│     [(0, 0), (1, 1)]
└     [(0, 0), (1, 1), (2, 1)]
[ Info: Interpolating the exponents in parameters..
┌ Info: Interpolating for degrees:
└ numerator 2, denominator 5
[ Info: Using 18 points..
[ Info: Using 36 points..
[ Info: Using 72 points..
[ Info: Success! 72 points used.
┌ Info: Output summary:
│     Maximal interpolated degrees are: 2 for num. and 5 for den.
│     Maximal number of interpolated terms are: 2 for num. and 3 for den.        
└ 
[ Info: Recovering the coefficients..
[ Info: Success! Used 1 prime in total :)
3-element Vector{AbstractAlgebra.Generic.MPoly{AbstractAlgebra.Generic.Frac{fmpq_mpoly}}}:
 y + (-a - 1)//(a^2*b^2 + a*b^4 + a*b^2)*z - 1//(a^2*b + a*b^3 + a*b)
 x + (-a - 1)//b*z
 z^2 + b//(a + 1)*z + b^2//(a + 1)
```
