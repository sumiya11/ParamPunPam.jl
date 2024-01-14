FIELDS = [Nemo.Native.GF(2^31 - 1), Nemo.Native.GF(2^62 + 135)]

@testset "Fast gcd" begin
    for ground in FIELDS
        R, x = polynomial_ring(ground, "x")
        cases = [
            (R(0), R(1)),
            (R(4), R(0)),
            (R(4), R(2)),
            (x, R(3)),
            (x, x^2),
            (x + 3, x + 2),
            (2x^3, 3x^4),
            (x + 1, 2(x + 1)),
            (x^2 + 1, x^2 + 2),
            (x^2 + 1, x^2 + 3),
            ((x - 1) * (x + 4) * (x - 8), (x - 8) * (x + 11)),
            ((x^4 + 4) * (x^2 + 3) * (x + 2), (x^2 + 2) * (x^4 + 4)),
            ((x - 2)^10 * (x^4 + 4)^5 * (x + 11), (x - 2)^7 * (x^2 + 2)^11 * (x + 11)),
            ((x - 1)^1000 * (x - 2)^200, (x - 1)^500 * (x - 2)^700)
        ]
        for (P, Q) in cases
            @test Nemo.gcd(P, Q) == ParamPunPam.fastgcd(P, Q)
        end

        # random tests
        for i in 1:100
            getrandpoly = () -> rand(R, 1:50, 1:5)
            a, b, c = getrandpoly(), getrandpoly(), getrandpoly()
            P = (1 + a * rand(0:1)) * (1 + b * rand(0:1)) * (1 + c * rand(0:1))
            Q = (1 + a * rand(0:1)) * (1 + b * rand(0:1)) * (1 + c * rand(0:1))
            @test Nemo.gcd(P, Q) == ParamPunPam.fastgcd(P, Q)
        end
    end
end

@testset "Pade approximation" begin
    for ground in FIELDS
        R, x = polynomial_ring(ground, "x")
        cases = [
            (R(1), R(0)),
            (R(4), R(0)),
            (R(4), R(2)),
            (x, R(3)),
            (x^5, x^5),
            (x^2, x),
            (2x^4, 3x^3),
            (x + 1, 2(x + 1)),
            (x^2 + 1, x^2 + 2),
            (x^2 + 1, x^2 + 3),
            ((x - 1) * (x + 4) * (x - 8), (x - 8) * (x + 11)),
            ((x^4 + 4) * (x^2 + 3) * (x + 2), (x^2 + 2) * (x^4 + 4)),
            ((x - 2)^20 * (x^4 + 4)^5 * (x + 11), (x - 2)^7 * (x^2 + 2)^11 * (x + 11)),
            ((x - 1)^1100 * (x - 2)^200, (x - 1)^500 * (x - 2)^700)
        ]
        for (P, Q) in cases
            dP, dQ = max(degree(P), 0), max(degree(Q), 0)
            dP == dQ && continue
            ks = rand(0:dP, 2)
            for k in ks
                r1, t1, s1 = ParamPunPam.constrainedEEA(P, Q, k)
                r2, t2, s2 = ParamPunPam.fastconstrainedEEA(P, Q, k)
                @test r1 == P * t1 + Q * s1
                @test r2 == P * t2 + Q * s2
                @test degree(r1) <= k
                @test degree(r2) <= k
                @test degree(r1) == degree(r2)
            end
        end

        # random tests
        for i in 1:100
            getrandpoly = () -> rand(R, 1:50, 1:5)
            a, b, c = getrandpoly(), getrandpoly(), getrandpoly()
            P = (1 + a * rand(0:1)) * (1 + b * rand(0:1)) * (1 + c * rand(0:1))
            Q = (1 + a * rand(0:1)) * (1 + b * rand(0:1)) * (1 + c * rand(0:1))
            dP, dQ = max(degree(P), 0), max(degree(Q), 0)
            if dP < dQ
                P, Q = Q, P
                dP, dQ = dQ, dP
            end
            dP == dQ && continue
            ks = rand(0:dP, 2)
            for k in ks
                r1, t1, s1 = ParamPunPam.constrainedEEA(P, Q, k)
                r2, t2, s2 = ParamPunPam.fastconstrainedEEA(P, Q, k)
                @test r1 == P * t1 + Q * s1
                @test r2 == P * t2 + Q * s2
                @test degree(r1) <= k
                @test degree(r2) <= k
                @test degree(r1) == degree(r2)
            end
        end
    end
end

@testset "Cauchy interpolation" begin
    for ground in FIELDS
        R, z = polynomial_ring(ground, "z")
        cases = [
            (R(0)) // (R(1)),
            (R(10)) // (R(1)),
            (R(-10)) // (R(1)),
            (z^2) // (z + 1),
            R(1) // (z + 8),
            R(-1) // (z + 8),
            (z - 1)^5 // (z + 1)^5,
            (z^3 - 3z^2 + 6z) // R(8),
            (z^3 - 3z^2 + 6z) // R(6),
            R(2) // R(3),
            ((z - 1)(z - 2)(z + 3)) // ((z + 6)^5),
            (3z - 1) // R(4),
            (z^8 + 8) // (z^4 + 2),
            (z^16 + 1) // (z^8 + 10),
            (z + 5)^10 // (z - 1)^7,
            (z)^100 // (z + 1)^100,
            R(19) // z^100,
            z^100 // R(15),
            5(z + 5)^2 // 17(z),
            (z + 18) // 18,
            (18z + 1) // 18,
            (-(z + 19) * (z + 17) + (z + 20) * (z + 21)) // R(1)
        ]
        evalfrac(f, x) = evaluate(numerator(f), x) // evaluate(denominator(f), x)
        for case in cases
            n, d = max(0, degree(numerator(case))), degree(denominator(case))
            c = ParamPunPam.CauchyInterpolator(R, n, d)
            xs = ParamPunPam.distinct_nonzero_points(ground, n + d + 2)
            ys = map(x -> evalfrac(case, x), xs)
            P, Q = ParamPunPam.interpolate!(c, xs, ys)
            @test isone(trailing_coefficient(Q))
            @test P // Q == case

            # Test for the case when n and d are upper bounds
            for _ in 1:5
                a, b = rand(1:3), rand(1:3)
                n, d = a * n + a, b * d + b
                subs = ground(rand(1:(2^20)))
                new_x = gen(R) * subs
                new_case = evaluate(case, new_x)
                c = ParamPunPam.CauchyInterpolator(R, n, d)
                xs = ParamPunPam.distinct_nonzero_points(ground, n + d + 2)
                ys = map(x -> evalfrac(new_case, x), xs)
                P, Q = ParamPunPam.interpolate!(c, xs, ys)
                @test isone(trailing_coefficient(Q))
                @test P // Q == new_case
            end
        end
    end
end
