# https://symbolicdata.github.io/PolynomialSystems
# https://github.com/symbolicdata/data

using Nemo

###
# Source: https://github.com/symbolicdata/data/blob/master/XMLResources/IntPS/Geometry.Chou.302_1.xml

P, (u1, u2, u3, u4, u5) = polynomial_ring(Nemo.QQ, [:u1, :u2, :u3, :u4, :u5])
R, (x1, x2, x3, x4, x5, x6, x7, x8) =
    polynomial_ring(fraction_field(P), [:x1, :x2, :x3, :x4, :x5, :x6, :x7, :x8])

system = [
    -x5 * u3 + x6 * u2,
    -x3 * u3 + x4 * u2,
    -x2 * u2 + x5 * u2 + x6 * u3 - u3 * u5,
    -x1 * u2 + x3 * u2 + x4 * u3 - u3 * u4,
    x2^2 * u3 - x2 * u1 * u3 + u1 * u2 * u5 - u2^2 * u5 - u3^2 * u5 + u3 * u5^2,
    x1 * x2 - x1 * x5 - x2 * x7 + x5 * x7 + x6 * x8 - x6 * u4,
    x1 * x2 - x2 * x3 - x1 * x7 + x3 * x7 + x4 * x8 - x4 * u5,
    x1^2 * u3 - x1 * u1 * u3 + u1 * u2 * u4 - u2^2 * u4 - u3^2 * u4 + u3 * u4^2
]

@time ParamPunPam.paramgb(system);

###
# source: https://github.com/symbolicdata/data/blob/master/XMLResources/IntPS/Geometry.Simson_3.xml

P, (u1, u2, u3, u4) = polynomial_ring(Nemo.QQ, [:u1, :u2, :u3, :u4])
R, (x1, x2, x3, x4, x5, x6, x7, x8, x9) = polynomial_ring(
    fraction_field(P),
    [:x1, :x2, :x3, :x4, :x5, :x6, :x7, :x8, :x9],
    ordering=:degrevlex
)

system = [
    x1 * x8 - x2 * x8 - x1 * u2 + x2 * u1 - x9 * u1 + x9 * u2,
    -x2 * x6 + x7 * u2,
    -x1 * x4 + x5 * u1,
    x3^2 + u3^2 - 2 * u3 * u4,
    -x2 * x3 + x2 * x7 + x6 * u2 - u2 * u3,
    -x1 * x3 + x1 * x5 + x4 * u1 - u1 * u3,
    x1 * x3 - x2 * x3 - x1 * x9 + x2 * x9 - x8 * u1 + x8 * u2 + u1 * u3 - u2 * u3,
    x2^2 + u2^2 - 2 * u2 * u4,
    x1^2 + u1^2 - 2 * u1 * u4
]

@time ParamPunPam.paramgb(system)

###
# source: https://github.com/symbolicdata/data/blob/master/XMLResources/IntPS/Vermeer.xml

using Groebner
R, (w, v, u, y, x) = polynomial_ring(QQ, [:w, :v, :u, :y, :x])
system = [
    v^2 + u^2 - 2 * v * y + y^2 - 2 * u * x + x^2 - 1,
    -u^3 + v^2,
    -3 * v * u^2 + 3 * u^2 * y - 2 * v * u + 2 * v * x,
    6 * w^2 * v * u^2 - 3 * w * u^2 - 2 * w * v + 1
]
@time groebner(system);
