using Singular

###
# Chou-302

P, (u1, u2, u3, u4, u5) = polynomial_ring(QQ, [:u1, :u2, :u3, :u4, :u5])
R, (x1, x2, x3, x4, x5, x6, x7, x8) =
    polynomial_ring(Singular.AbstractAlgebra.fraction_field(P), [:x1, :x2, :x3, :x4, :x5, :x6, :x7, :x8], ordering=:degrevlex)

system = [
    -x5 * (u3//1) + x6 * (u2//1),
    -x3 * (u3//1) + x4 * (u2//1),
    -x2 * (u2//1) + x5 * (u2//1) + x6 * (u3//1) - u3 * (u5//1),
    -x1 * (u2//1) + x3 * (u2//1) + x4 * (u3//1) - u3 * (u4//1),
    x2^2 * (u3//1) - x2 * (u1//1) * (u3//1) + u1 * (u2//1) * (u5//1) - u2^2 * (u5//1) - u3^2 * (u5//1) + u3 * (u5//1)^2,
    x1 * x2 - x1 * x5 - x2 * x7 + x5 * x7 + x6 * x8 - x6 * (u4//1),
    x1 * x2 - x2 * x3 - x1 * x7 + x3 * x7 + x4 * x8 - x4 * (u5//1),
    x1^2 * (u3//1) - x1 * (u1//1) * (u3//1) + u1 * (u2//1) * (u4//1) - u2^2 * (u4//1) - u3^2 * (u4//1) + u3 * (u4//1)^2
]

@time Singular.slimgb(Singular.Ideal(R, system));
# 54.250032 seconds (12.79 M allocations: 33.071 GiB, 0.57% gc time, 0.04% compilation time)

###
# Simson-3

P, (u1, u2, u3, u4) = polynomial_ring(QQ, [:u1, :u2, :u3, :u4])
R, (x1, x2, x3, x4, x5, x6, x7, x8, x9) = polynomial_ring(
    Singular.AbstractAlgebra.fraction_field(P),
    [:x1, :x2, :x3, :x4, :x5, :x6, :x7, :x8, :x9],
    ordering=:degrevlex
)

system = [
    x1 * x8 - x2 * x8 - x1 * (u2//1) + x2 * (u1//1) - x9 * (u1//1) + x9 * (u2//1),
    -x2 * x6 + x7 * (u2//1),
    -x1 * x4 + x5 * (u1//1),
    x3^2 + (u3^2//1) - 2 * (u3//1) * (u4//1),
    -x2 * x3 + x2 * x7 + x6 * (u2//1) - u2 * (u3//1),
    -x1 * x3 + x1 * x5 + x4 * (u1//1) - u1 * (u3//1),
    x1 * x3 - x2 * x3 - x1 * x9 + x2 * x9 - x8 * (u1//1) + x8 * (u2//1) + u1 * (u3//1) - u2 * (u3//1),
    x2^2 + (u2^2//1) - 2 * (u2//1) * (u4//1),
    x1^2 + (u1^2//1) - 2 * (u1//1) * (u4//1)
]

@time Singular.slimgb(Singular.Ideal(R, system));
# Timeout (>1 minute)

###
# Param-1 (katsura)

P, (a,b) = polynomial_ring(QQ, [:a,:b])
R, (x0,x1,x2,x3,x4) = polynomial_ring(Singular.AbstractAlgebra.fraction_field(P), [:x0,:x1,:x2,:x3,:x4], ordering=:degrevlex)
system = [
    x0^2 + 2*x1^2 + 2*x2^2 + 2*x3^2 + 2*x4^2 - x0 - b//1,
    2*x0*x1 + 2*x1*x2 + 2*x2*x3 + 2*x3*x4 - x1,
    x1^2 + 2*x0*x2 + 2*x1*x3 + 2*x2*x4 - x2,
    2*x1*x2 + 2*x0*x3 + 2*x1*x4 - x3,
    b//1*x0 + 2*x1 + 2*x2 + 2*x3 + 2*x4 - a//1
]

@time Singular.slimgb(Singular.Ideal(R, system));
# Timeout (>1 minute)

###
# Param-2 (katsura)

system = [
    (b+1)//1*(x0^2 + 2*x1^2 + 2*x2^2 + 2*x3^2 + 2*x4^2 - x0) - b//1,
    (a+b)//1*(2*x0*x1 + 2*x1*x2 + 2*x2*x3 + 2*x3*x4 - x1) - b//1,
    b//1*(x1^2 + 2*x0*x2 + 2*x1*x3 + 2*x2*x4 - x2) - a//1,
    (a-1)//1*(2*x1*x2 + 2*x0*x3 + 2*x1*x4 - x3) - b//1,
    x0 + 2*x1 + 2*x2 + 2*x3 + 2*x4 - a^3//1
]

@time Singular.slimgb(Singular.Ideal(R, system));
# Timeout (>1 minute)

nothing

