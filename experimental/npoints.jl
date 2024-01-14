using Nemo

Rparam, (a1, a2, a3, a4, a5, a6) = polynomial_ring(Nemo.QQ, ["a$i" for i in 1:6])
R, (x1, x2, x3) =
    polynomial_ring(Nemo.fraction_field(Rparam), ["x1", "x2", "x3"], ordering=:degrevlex)
f = [
    x1 * x2^2 + (a1 + a2 + a3) * x1 * x3^2 - (a2 + a5 + a6),
    a3 * a5 * a6^2 * x1^2 * x2 + (a1 + a2) * x2 * x3^2 - a2 * x2 + a1,
    (a1 + 3a2) * x1^2 * x3 + a1 * x2^2 * x3 - a2 * x3 + (a1 + a3 + a6)
]

@profview for _ in 1:100
    ParamPunPam.paramgb(f, up_to_degree=(2, 2))
end

f = [x1 * x2^2 + (a1^3 + a2 + a3) // (a1 + a3^9 + a6)]
ParamPunPam.paramgb(f, up_to_degree=(2, 2))

ParamPunPam._runtime_data
