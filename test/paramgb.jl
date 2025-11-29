using Random, Logging

Groebner = ParamPunPam.Groebner
interpolators_to_test = [:VanDerHoevenLecerf, :CuytLee]

function test_paramgb(cases, answers; kwargs...)
    for (case, answer) in zip(cases, answers)
        gb = ParamPunPam.paramgb(case; kwargs...)
        # Test that coefficients up to the given degrees coincide 
        if !haskey(kwargs, :up_to_degree)
            @test gb == answer
        else
            up_to_degree = get(kwargs, :up_to_degree, (0, 0))
            # `get` always succeeds
            @assert !any(iszero.(up_to_degree))
            @test length(gb) == length(answer)
            for i in 1:length(gb)
                @test length(gb[i]) == length(answer[i])
                for j in 1:length(gb[i])
                    true_c = coeff(answer[i], j)
                    gb_c = coeff(gb[i], j)
                    dn, dd = map(total_degree, (numerator(true_c), denominator(true_c)))
                    if all((dn, dd) .<= up_to_degree)
                        @test true_c == gb_c
                    else
                        @test isone(gb_c)
                    end
                end
            end
        end
    end
end

@testset "GB over Q(a...)" begin
    Ra, (a, b) = polynomial_ring(Nemo.QQ, ["a", "b"])
    Rx, (x, y, z) = polynomial_ring(Nemo.fraction_field(Ra), ["x", "y", "z"])
    @test ParamPunPam.paramgb([Rx(2)]) == [Rx(1)]
    @test_broken ParamPunPam.paramgb([Rx(0)]) == [Rx(1)]

    for interpolator in interpolators_to_test
        for param_ord in [:lex, :deglex, :degrevlex]
            for up_to_degree in [(Inf, Inf), (1, 1), (2, 2), (4, 4), (8, 8), (16, 16)]
                Ra, (a, b) = polynomial_ring(Nemo.QQ, ["a", "b"], internal_ordering=param_ord)
                Rx, (x, y, z) = polynomial_ring(Nemo.fraction_field(Ra), ["x", "y", "z"], internal_ordering=:degrevlex)
                # test that invalid keyword arguments are reported.
                # Notice `up_to_degreeS`.
                @test_throws AssertionError ParamPunPam.paramgb([x], up_to_degrees=(3, 3)) == [x]

                cases = [
                    [x, x + a],
                    [(x + a)^2],
                    [x + a^2, x * y + a^10],
                    [x^2 + a^3 + 1],
                    [x^2 + a^2 + b^2 + a * b + 9],
                    [x^2 - x + a * y^2 + a * z^2, a * x * y + a * y * z - y, x + a * y + a * z - 1]
                ]
                answers = [
                    [Rx(1)],
                    [(x + a)^2],
                    [y - a^8, x + a^2],
                    [x^2 + a^3 + 1],
                    [x^2 + a^2 + b^2 + a * b + 9],
                    [
                        x + a * y + a * z - 1,
                        y * z + (a^2 + a) // (a^2 + 1) * z^2 - 1 // (a^3 + a) * y - a // (a^2 + 1) * z,
                        y^2 + (-a^2 + 1) // (a^2 + 1) * z^2 + (-a + 1) // (a^2 + 1) * y + (a - 1) // (a^2 + 1) * z,
                        z^3 +
                        (-3 // 2 * a^5 + a^4 - a^3 + 1 // 2 * a^2 - 1 // 2 * a - 1 // 2) //
                        (a^6 + 1 // 2 * a^5 + a^4 + a^3 + 1 // 2 * a) * z^2 +
                        (1 // 2 * a^3 - 1 // 2) // (a^6 + 1 // 2 * a^5 + a^4 + a^3 + 1 // 2 * a) * y +
                        (1 // 2 * a^4 - a^3 + 1 // 2 * a^2 - 1 // 2 * a + 1 // 2) //
                        (a^6 + 1 // 2 * a^5 + a^4 + a^3 + 1 // 2 * a) * z
                    ]
                ]
                test_paramgb(cases, answers, rational_interpolator=interpolator)
                test_paramgb(cases, answers, up_to_degree=up_to_degree, rational_interpolator=interpolator)

                Ra, (a1, a2, a3, a4, a5) =
                    polynomial_ring(Nemo.QQ, ["a1", "a2", "a3", "a4", "a5"], internal_ordering=param_ord)
                a = [a1, a2, a3, a4, a5]
                Rx, (x, y, z) = polynomial_ring(Nemo.fraction_field(Ra), ["x", "y", "z"], internal_ordering=:degrevlex)
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
                test_paramgb(cases, answers, rational_interpolator=interpolator)
                test_paramgb(cases, answers, up_to_degree=up_to_degree, rational_interpolator=interpolator)

                Ra, (a, b, c) = polynomial_ring(Nemo.QQ, ["a", "b", "c"], internal_ordering=param_ord)
                Rx, (x, y, z) = polynomial_ring(Nemo.fraction_field(Ra), ["x", "y", "z"], internal_ordering=:deglex)
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
                test_paramgb(cases, answers, rational_interpolator=interpolator)
                test_paramgb(cases, answers, up_to_degree=up_to_degree, rational_interpolator=interpolator)
            end
        end
    end
end

@testset "Monomial orderings" begin
    for interpolator in interpolators_to_test
        Rparam, (a, b, c) = polynomial_ring(Nemo.QQ, ["a", "b", "c"])
        R, (x1, x2, x3) = polynomial_ring(Nemo.fraction_field(Rparam), ["x1", "x2", "x3"], internal_ordering=:degrevlex)

        f = [a * x1 - b * x2 + c * x3 - 1]
        gb1 = ParamPunPam.paramgb(f, ordering=ParamPunPam.Lex())
        @test gb1 == [x1 - (b // a) * x2 + (c // a) * x3 - 1 // a]
        gb2 = ParamPunPam.paramgb(f, ordering=ParamPunPam.Lex(x2, x1, x3))
        @test gb2 == [-(a // b) * x1 + x2 - (c // b) * x3 + 1 // b]

        f = [x1 - a * b^2 + a^2 * b + a^3 + b^3, x3 - (a * b^3 + b * a^3)^3, x2 - a * b^3 + a^2 * b + a^4 + b^4]
        gb = ParamPunPam.paramgb(f, ordering=ParamPunPam.Lex(x3, x2, x1))
        @test gb == [x1 - a * b^2 + a^2 * b + a^3 + b^3, x2 - a * b^3 + a^2 * b + a^4 + b^4, x3 - (a * b^3 + b * a^3)^3]
	degs = paramgb_only_degrees(f, ordering=ParamPunPam.Lex(x3, x2, x1), up_to_degree=(5, 10))
	@test degs == [[(0, 0), (3, 0)], [(0, 0), (4, 0)], [(0, 0), (-1, -1)]]

        f = [x2 + a, x1 + b, x3 + c]
        gb1 = ParamPunPam.paramgb(f, ordering=ParamPunPam.Lex())
        gb2 = ParamPunPam.paramgb(f, ordering=ParamPunPam.DegLex())
        gb3 = ParamPunPam.paramgb(f, ordering=ParamPunPam.DegRevLex())
        @test gb1 == gb2 == gb3 == ParamPunPam.paramgb(f)
        @test parent(gb1[1]) == R

        gb1 = ParamPunPam.paramgb(f, ordering=ParamPunPam.Lex(x3, x2, x1))
        @test gb1 == [x1 + b, x2 + a, x3 + c]

        gb1 = ParamPunPam.paramgb(f, ordering=ParamPunPam.Lex(x2, x1, x3))
        @test gb1 == [x3 + c, x1 + b, x2 + a]

        gb1 = ParamPunPam.paramgb(f, ordering=ParamPunPam.DegRevLex(x3, x2, x1))
        @test gb1 == [x1 + b, x2 + a, x3 + c]

        # The order of output persists
        for ord in [:lex, :deglex, :degrevlex]
            Rparam, (a, b, c) = polynomial_ring(Nemo.QQ, ["a", "b", "c"])
            R, (x1, x2, x3) = polynomial_ring(Nemo.fraction_field(Rparam), ["x1", "x2", "x3"], internal_ordering=ord)

            cases = ([x1, x2, x3], [x1 + 2a, c * x2 + 3b, x3], [x1 + x2 + x3, x1 + x2, x1])
            for (gb_ord, var_to_index) in [
                (ParamPunPam.Lex(), Dict(x1 => 3, x2 => 2, x3 => 1)),
                (ParamPunPam.DegLex(), Dict(x1 => 3, x2 => 2, x3 => 1)),
                (ParamPunPam.DegRevLex(), Dict(x1 => 3, x2 => 2, x3 => 1)),
                (ParamPunPam.Lex(x2, x1, x3), Dict(x1 => 2, x2 => 3, x3 => 1)),
                (ParamPunPam.Lex(x3, x2, x1), Dict(x1 => 1, x2 => 2, x3 => 3))
            ]
                for case in cases
                    gb = ParamPunPam.paramgb(case, ordering=gb_ord)
                    for (rk, f) in enumerate(gb)
                        m = Nemo.leading_monomial(f)
                        @test rk == var_to_index[m]
                    end
                end
            end
        end

        # GB with no parameters coincides with the numerical GB
        ord = :degrevlex
        cases = [
            Groebner.Examples.noonn(3, internal_ordering=ord),
            Groebner.Examples.noonn(4, internal_ordering=ord),
            Groebner.Examples.katsuran(3, internal_ordering=ord),
            Groebner.Examples.katsuran(3, internal_ordering=ord),
            Groebner.Examples.cyclicn(3, internal_ordering=ord)
        ]
        for case in cases
            Rparam, (a, b, c) = polynomial_ring(Nemo.QQ, ["a", "b", "c"])
            xs = gens(parent(case[1]))
            (R, xs_frac) = polynomial_ring(Nemo.fraction_field(Rparam), map(repr, xs), internal_ordering=ord)
            gb_ords = []
            append!(gb_ords, [Groebner.DegLex(), Groebner.DegLex()])
            append!(gb_ords, [Groebner.DegLex(Random.shuffle(xs)) for _ in 1:10])
            append!(gb_ords, [Groebner.DegRevLex(Random.shuffle(xs)) for _ in 1:10])
            for gb_ord in gb_ords
                case_frac = map(f -> evaluate(f, xs_frac), case)

                gb_numeric = Groebner.groebner(case, ordering=gb_ord)
                gb_parametric = ParamPunPam.paramgb(case_frac, ordering=gb_ord, rational_interpolator=interpolator)

                gb_numeric_frac = map(f -> evaluate(f, xs_frac), gb_numeric)

                @test gb_numeric_frac == gb_parametric
            end
        end
    end
end

@testset "Multi-modular over the rationals" begin
    Rparam, (a, b, c) = polynomial_ring(Nemo.QQ, ["a", "b", "c"])
    R, (x1, x2, x3) = polynomial_ring(Nemo.fraction_field(Rparam), ["x1", "x2", "x3"], internal_ordering=:degrevlex)
    cases = [
        [x2^2 + (3 // 23) * x2 + (4 // 25) * x3, x1 + 5 // 7],
        [x1 + BigInt(2)^100, x2^2 - BigInt(2^31 + 1)^5],
        [x1 + (2^20 * a + c) * x3 - 2^21 * b^2, x2 + 2^19],
        [x1 + (2^40 * a + c) * x3 - 2^41 * b^2, x2 + 2^39],
        [x1 + (BigInt(2)^80 * a + c) * x3 - BigInt(2)^81 * b^2, x2 + BigInt(2)^79],
        [x2^2 - x2 + 1, (2^30) * x1 + (2^31 + 5) * x2],
        [x1 + BigInt(2)^1000]
    ]
    answers = [
        [x1 + 5 // 7, x2^2 + (3 // 23) * x2 + (4 // 25) * x3],
        [x1 + BigInt(2)^100, x2^2 - BigInt(2^31 + 1)^5],
        [x2 + 2^19, x1 + (2^20 * a + c) * x3 - 2^21 * b^2],
        [x2 + 2^39, x1 + (2^40 * a + c) * x3 - 2^41 * b^2],
        [x2 + BigInt(2)^79, x1 + (BigInt(2)^80 * a + c) * x3 - BigInt(2)^81 * b^2],
        [x1 + (2^31 + 5) // (2^30) * x2, x2^2 - x2 + 1],
        [x1 + BigInt(2)^1000]
    ]
    test_paramgb(cases, answers)

    F = [x1 + a * b^3 // BigInt(2)^100 * x2]
    ParamPunPam.paramgb(F)
end

@testset "Noon" begin
    Rparam, (a, b, c) = polynomial_ring(Nemo.QQ, ["a", "b", "c"])
    R, (x1, x2, x3) = polynomial_ring(Nemo.fraction_field(Rparam), ["x1", "x2", "x3"], internal_ordering=:degrevlex)
    f = [
        a * x1 * x2^2 + (a + c) * x1 * x3^2 - b * x1 + a,
        a * x1^2 * x2 + (a + b) * x2 * x3^2 - b * x2 + a,
        a * x1^2 * x3 + a * x2^2 * x3 - b * x3 + a
    ]
    ParamPunPam.paramgb(f)

    Rparam, (a1, a2, a3, a4, a5, a6) = polynomial_ring(Nemo.QQ, ["a$i" for i in 1:6])
    R, (x1, x2, x3) = polynomial_ring(Nemo.fraction_field(Rparam), ["x1", "x2", "x3"], internal_ordering=:degrevlex)
    f = [
        x1 * x2^2 + (a1 + a2 + a3) * x1 * x3^2 - (a2 + a5 + a6),
        x1^2 * x2 + (a1 + a2) * x2 * x3^2 - a2 * x2 + a1,
        x1^2 * x3 + a1 * x2^2 * x3 - a2 * x3 + (a1 + a3 + a6)
    ]
    ParamPunPam.paramgb(f)

    Rparam, Ai = polynomial_ring(Nemo.QQ, ["A$i" for i in 1:7])
    R, (x, y, z) = polynomial_ring(Nemo.fraction_field(Rparam), ["x", "y", "z"], internal_ordering=:degrevlex)
    f = [sum(Ai) * x * y + (Ai[1] + Ai[6]) // (sum(Ai)), y * z - sum(Ai), x^2 + Ai[5] // (sum(Ai))]
    ParamPunPam.paramgb(f)
end
