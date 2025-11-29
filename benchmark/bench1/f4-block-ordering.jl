using Nemo, Groebner

###
# Chou-302

R, (u1, u2, u3, u4, u5, x1, x2, x3, x4, x5, x6, x7, x8) =
    polynomial_ring(QQ, [:u1, :u2, :u3, :u4, :u5, :x1, :x2, :x3, :x4, :x5, :x6, :x7, :x8])

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

@time gb = groebner(system, ordering=DegRevLex(x1, x2, x3, x4, x5, x6, x7, x8)*DegRevLex(u1, u2, u3, u4, u5));
#   2.646309 seconds (5.24 M allocations: 495.906 MiB, 4.14% gc time)

###
# Simson-3

R, (u1, u2, u3, u4, x1, x2, x3, x4, x5, x6, x7, x8, x9) = polynomial_ring(
    QQ,
    [:u1, :u2, :u3, :u4, :x1, :x2, :x3, :x4, :x5, :x6, :x7, :x8, :x9],
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

@time gb = Groebner.groebner(system, ordering=DegRevLex(x1, x2, x3, x4, x5, x6, x7, x8, x9)*DegRevLex(u1, u2, u3, u4))
#   1.865856 seconds (3.80 M allocations: 345.477 MiB, 6.97% gc time)

###
# Param-1 (katsura)

R, (a,b,x0,x1,x2,x3,x4) = polynomial_ring(QQ, [:a,:b,:x0,:x1,:x2,:x3,:x4], internal_ordering=:degrevlex)
system = [
    x0^2 + 2*x1^2 + 2*x2^2 + 2*x3^2 + 2*x4^2 - x0 - b,
    2*x0*x1 + 2*x1*x2 + 2*x2*x3 + 2*x3*x4 - x1,
    x1^2 + 2*x0*x2 + 2*x1*x3 + 2*x2*x4 - x2,
    2*x1*x2 + 2*x0*x3 + 2*x1*x4 - x3,
    b*x0 + 2*x1 + 2*x2 + 2*x3 + 2*x4 - a
]

@time gb = groebner(system, ordering=DegRevLex(x0,x1,x2,x3,x4)*DegRevLex(a,b));
#   2.604742 seconds (4.26 M allocations: 522.964 MiB, 4.05% gc time)

###
# Param-2 (katsura)

system = [
    (b+1)*(x0^2 + 2*x1^2 + 2*x2^2 + 2*x3^2 + 2*x4^2 - x0) - b,
    (a+b)*(2*x0*x1 + 2*x1*x2 + 2*x2*x3 + 2*x3*x4 - x1) - b,
    b*(x1^2 + 2*x0*x2 + 2*x1*x3 + 2*x2*x4 - x2) - a,
    (a-1)*(2*x1*x2 + 2*x0*x3 + 2*x1*x4 - x3) - b,
    x0 + 2*x1 + 2*x2 + 2*x3 + 2*x4 - a^3
]

@time gb = groebner(system, ordering=DegRevLex(x0,x1,x2,x3,x4)*DegRevLex(a,b));
#  9.152303 seconds (9.28 M allocations: 1.570 GiB, 3.07% gc time)

nothing
