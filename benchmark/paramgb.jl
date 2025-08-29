using ParamPunPam, Nemo

###
# Chou-302

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

@time gb = ParamPunPam.paramgb(system, ordering=DegRevLex());
#  1.981330 seconds (12.26 M allocations: 660.422 MiB, 17.96% gc time, 0.63% compilation time)

###
# Simson-3

P, (u1, u2, u3, u4) = polynomial_ring(Nemo.QQ, [:u1, :u2, :u3, :u4])
R, (x1, x2, x3, x4, x5, x6, x7, x8, x9) = polynomial_ring(
    fraction_field(P),
    [:x1, :x2, :x3, :x4, :x5, :x6, :x7, :x8, :x9],
    internal_ordering=:degrevlex
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

@time gb = ParamPunPam.paramgb(system, ordering=DegRevLex());
#   0.788609 seconds (2.94 M allocations: 167.057 MiB, 25.87% gc time, 2.83% compilation time)

###
# Param-1 (katsura)

P, (a,b) = polynomial_ring(QQ, ["a","b"])
R, (x0,x1,x2,x3,x4) = polynomial_ring(fraction_field(P), [:x0,:x1,:x2,:x3,:x4], internal_ordering=:degrevlex)
system = [
    x0^2 + 2*x1^2 + 2*x2^2 + 2*x3^2 + 2*x4^2 - x0 - b,
    2*x0*x1 + 2*x1*x2 + 2*x2*x3 + 2*x3*x4 - x1,
    x1^2 + 2*x0*x2 + 2*x1*x3 + 2*x2*x4 - x2,
    2*x1*x2 + 2*x0*x3 + 2*x1*x4 - x3,
    b*x0 + 2*x1 + 2*x2 + 2*x3 + 2*x4 - a
]

@time gb = ParamPunPam.paramgb(system, ordering=DegRevLex());
# Timeout (>1 minute)

###
# Param-2 (katsura)

system = [
    (b+1)*(x0^2 + 2*x1^2 + 2*x2^2 + 2*x3^2 + 2*x4^2 - x0) - b,
    (a+b)*(2*x0*x1 + 2*x1*x2 + 2*x2*x3 + 2*x3*x4 - x1) - b,
    b*(x1^2 + 2*x0*x2 + 2*x1*x3 + 2*x2*x4 - x2) - a,
    (a-1)*(2*x1*x2 + 2*x0*x3 + 2*x1*x4 - x3) - b,
    x0 + 2*x1 + 2*x2 + 2*x3 + 2*x4 - a^3
]

@time gb = ParamPunPam.paramgb(system, ordering=DegRevLex());
# 13.389742 seconds (63.77 M allocations: 4.118 GiB, 12.63% gc time, 0.89% compilation time)

nothing
