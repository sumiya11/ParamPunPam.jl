
function test_paramgb(cases, answers; kwargs...)
    for (case, answer) in zip(cases, answers)
        gb = ParamPunPam.paramgb(case; kwargs...)
        @info "" case gb kwargs
        # Test that coefficients up to the given degrees coincide 
        if haskey(kwargs, :up_to_degree)
            up_to_degree = get(kwargs, :up_to_degree, (0, 0))
            @assert !any(iszero.(up_to_degree))
            @test length(gb) == length(answer)
            for i in 1:length(gb)
                @test length(gb[i]) == length(answer[i])
                for j in 1:length(gb[i])
                    true_c = coeff(answer[i], j)
                    gb_c = coeff(gb[i], j)
                    dn, dd = map(total_degree, (numerator(true_c), denominator(true_c)))
                    if all((dn, dd) .< up_to_degree)
                        @test true_c == gb_c
                    end
                end
            end
        else
            @test gb == answer
        end
    end
end

@testset "GB over Q(a...)" begin
    for param_ord in [:lex, :deglex, :degrevlex]
        for up_to_degree in [(Inf, Inf), (1, 1), (2, 2), (4, 4), (8, 8), (16, 16)]
            Ra, (a,) = PolynomialRing(Nemo.QQ, ["a"], ordering=param_ord)
            Rx, (x, y, z) =
                PolynomialRing(Nemo.FractionField(Ra), ["x", "y", "z"], ordering=:degrevlex)
            cases = [
                [x, x + a],
                [(x + a)^2],
                [x + a^2, x * y + a^10],
                [
                    x^2 - x + a * y^2 + a * z^2,
                    a * x * y + a * y * z - y,
                    x + a * y + a * z - 1
                ]
            ]
            answers = [
                [Rx(1)],
                [(x + a)^2],
                [y - a^8, x + a^2],
                [
                    x + a * y + a * z - 1,
                    y * z + (a^2 + a) // (a^2 + 1) * z^2 - 1 // (a^3 + a) * y -
                    a // (a^2 + 1) * z,
                    y^2 +
                    (-a^2 + 1) // (a^2 + 1) * z^2 +
                    (-a + 1) // (a^2 + 1) * y +
                    (a - 1) // (a^2 + 1) * z,
                    z^3 +
                    (-3 // 2 * a^5 + a^4 - a^3 + 1 // 2 * a^2 - 1 // 2 * a - 1 // 2) //
                    (a^6 + 1 // 2 * a^5 + a^4 + a^3 + 1 // 2 * a) * z^2 +
                    (1 // 2 * a^3 - 1 // 2) //
                    (a^6 + 1 // 2 * a^5 + a^4 + a^3 + 1 // 2 * a) * y +
                    (1 // 2 * a^4 - a^3 + 1 // 2 * a^2 - 1 // 2 * a + 1 // 2) //
                    (a^6 + 1 // 2 * a^5 + a^4 + a^3 + 1 // 2 * a) * z
                ]
            ]
            test_paramgb(cases, answers)
            test_paramgb(cases, answers, up_to_degree=up_to_degree)

            Ra, (a1, a2, a3, a4, a5) =
                PolynomialRing(Nemo.QQ, ["a1", "a2", "a3", "a4", "a5"], ordering=param_ord)
            a = [a1, a2, a3, a4, a5]
            Rx, (x, y, z) =
                PolynomialRing(Nemo.FractionField(Ra), ["x", "y", "z"], ordering=:degrevlex)
            cases = [
                [x, y + 1, z + 2],
                [x + a1 * a2 - a3 * a4],
                [x + a1, x + a2, x + a3],
                [x + a1 // a2, y - (a1 - 1) // (a2 - 1), z + a4 // a3 * a5],
                [x + prod(a) + sum(a) - 1],
                [(x - a1) * (y - a2) * (z - a3) * (x - a4) * (x - a5)]
            ]
            answers = [
                [z + 2, y + 1, x],
                [x + a1 * a2 - a3 * a4],
                [Rx(1)],
                [z + a4 // a3 * a5, y - (a1 - 1) // (a2 - 1), x + a1 // a2],
                [x + prod(a) + sum(a) - 1],
                [(x - a1) * (y - a2) * (z - a3) * (x - a4) * (x - a5)]
            ]
            test_paramgb(cases, answers)
            test_paramgb(cases, answers, up_to_degree=up_to_degree)

            Ra, (a, b, c) = PolynomialRing(Nemo.QQ, ["a", "b", "c"], ordering=param_ord)
            Rx, (x, y, z) =
                PolynomialRing(Nemo.FractionField(Ra), ["x", "y", "z"], ordering=:deglex)
            cases = [
                [x + (a + b + c)^3 // (a * c * b)^2],
                [x + a // b, y + b // a, z + (a + b) // c],
                [z + (a + b + c)^8 // c, y + b // a, x + a // b],
                [z + a^25 // c, y + b^25 // a, x + a // b]
            ]
            answers = [
                [x + (a + b + c)^3 // (a * c * b)^2],
                [z + (a + b) // c, y + b // a, x + a // b],
                [z + (a + b + c)^8 // c, y + b // a, x + a // b],
                [z + a^25 // c, y + b^25 // a, x + a // b]
            ]
            test_paramgb(cases, answers)
            test_paramgb(cases, answers, up_to_degree=up_to_degree)
        end
    end
end

@testset "Multi-modular" begin
    Rparam, (a, b, c) = PolynomialRing(Nemo.QQ, ["a", "b", "c"])
    R, (x1, x2, x3) =
        PolynomialRing(Nemo.FractionField(Rparam), ["x1", "x2", "x3"], ordering=:degrevlex)
    cases = [
        [x2^2 + (3 // 23) * x2 + (4 // 25) * x3, x1 + 5 // 7],
        [x1 + BigInt(2)^100, x2^2 - BigInt(2^31 + 1)^5],
        [x1 + (2^20 * a + c) * x3 - 2^21 * b^2, x2 + 2^19],
        [x1 + (2^40 * a + c) * x3 - 2^41 * b^2, x2 + 2^39],
        [x1 + (BigInt(2)^80 * a + c) * x3 - BigInt(2)^81 * b^2, x2 + BigInt(2)^79],
        [x2^2 - x2 + 1, (2^30) * x1 + (2^31 + 5) * x2]
    ]
    answers = [
        [x1 + 5 // 7, x2^2 + (3 // 23) * x2 + (4 // 25) * x3],
        [x1 + BigInt(2)^100, x2^2 - BigInt(2^31 + 1)^5],
        [x2 + 2^19, x1 + (2^20 * a + c) * x3 - 2^21 * b^2],
        [x2 + 2^39, x1 + (2^40 * a + c) * x3 - 2^41 * b^2],
        [x2 + BigInt(2)^79, x1 + (BigInt(2)^80 * a + c) * x3 - BigInt(2)^81 * b^2],
        [(2^30) * x1 + (2^31 + 5) * x2, x2^2 - x2 + 1]
    ]
    # test_paramgb(cases, answers)
    @test_broken false
end

@testset "Noon" begin
    Rparam, (a, b, c) = PolynomialRing(Nemo.QQ, ["a", "b", "c"])
    R, (x1, x2, x3) =
        PolynomialRing(Nemo.FractionField(Rparam), ["x1", "x2", "x3"], ordering=:degrevlex)
    f = [
        a * x1 * x2^2 + (a + c) * x1 * x3^2 - b * x1 + a,
        a * x1^2 * x2 + (a + b) * x2 * x3^2 - b * x2 + a,
        a * x1^2 * x3 + a * x2^2 * x3 - b * x3 + a
    ]
    ParamPunPam.paramgb(f)

    Rparam, (a1, a2, a3, a4, a5, a6) = PolynomialRing(Nemo.QQ, ["a$i" for i in 1:6])
    R, (x1, x2, x3) =
        PolynomialRing(Nemo.FractionField(Rparam), ["x1", "x2", "x3"], ordering=:degrevlex)
    f = [
        x1 * x2^2 + (a1 + a2 + a3) * x1 * x3^2 - (a2 + a5 + a6),
        x1^2 * x2 + (a1 + a2) * x2 * x3^2 - a2 * x2 + a1,
        x1^2 * x3 + a1 * x2^2 * x3 - a2 * x3 + (a1 + a3 + a6)
    ]
    ParamPunPam.paramgb(f)

    Rparam, Ai = PolynomialRing(Nemo.QQ, ["A$i" for i in 1:7])
    R, (x, y, z) =
        PolynomialRing(Nemo.FractionField(Rparam), ["x", "y", "z"], ordering=:degrevlex)
    f = [
        sum(Ai) * x * y + (Ai[1] + Ai[6]) // (sum(Ai)),
        y * z - sum(Ai),
        x^2 + Ai[5] // (sum(Ai))
    ]
    ParamPunPam.paramgb(f)
end
