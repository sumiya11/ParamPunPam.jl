## Benchmarks

29 August, 2025. Alexander Demin.

Benchmarking the computation of Groebner bases in $\mathbb{Q}(x_1,\ldots,x_m)[y_1,\ldots,y_n]$.

#### Benchmark scripts

> f4-block-ordering.jl

Method: Groebner.jl, computing in $\mathbb{Q}[x_1,\ldots,x_m, y_1,\ldots,y_n]$ using a block ordering with $x_1,\ldots,x_m < y_1,\ldots,y_n$.

> f4-direct.jl

Method: Groebner.jl, directly computing in $\mathbb{Q}(x_1,\ldots,x_m)[y_1,\ldots,y_n]$.

> paramgb.jl

Method: ParamPunPam.jl, computing in $\mathbb{Q}(x_1,\ldots,x_m)[y_1,\ldots,y_n]$ using sparse interpolation.

> slimgb.jl

Method: slimgb from Singular.jl.

#### Benchmark systems
- [Simson-3](https://github.com/symbolicdata/data/blob/master/XMLResources/IntPS/Geometry.Simson_3.xml)
- [Chou-302](https://github.com/symbolicdata/data/blob/master/XMLResources/IntPS/Geometry.Chou.302_1.xml)
- Param-1
- Param-2
- Goodwin (the source is not available)

#### References
- https://github.com/symbolicdata/data
- https://symbolicdata.github.io/PolynomialSystems
