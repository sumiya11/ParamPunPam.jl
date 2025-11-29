FIELDS = [Nemo.Native.GF(2^62 + 135), Nemo.Native.GF(2^31 - 1)]

t1,t2 = 0, 0

@testset "Univariate interpolate" begin
    for ground in FIELDS
        R, x = polynomial_ring(ground, "x")
        cases = [
            R(0),
            R(4),
            x,
            x^2,
            2x^3,
            3x^4,
            x + 2,
            x^2 + 1,
            (x - 1) * (x + 4) * (x - 8),
            (x^4 + 4) * (x^2 + 3) * (x + 2),
            (x - 2)^10 * (x^4 + 4)^5 * (x + 11),
            (x - 1)^10 * (x - 2)^20,
            (x - 1)^50 * (x - 2)^70,
            x^20 + 1,
            10(x + 20)^100
        ]
        for f in cases
            xs = ParamPunPam.distinct_nonzero_points(ground, max(degree(f) + 1, 1))
            ys = map(x -> evaluate(f, x), xs)
            global t1 += @elapsed Nemo.interpolate(R, xs, ys)
            global t2 += @elapsed ParamPunPam.fastpolyinterpolate(R, xs, ys)
            @test f == Nemo.interpolate(R, xs, ys) == ParamPunPam.fastpolyinterpolate(R, xs, ys)
        end

        # random tests
        for i in 1:100
            getrandpoly = () -> rand(R, 1:50, 1:5)
            a, b, c = getrandpoly(), getrandpoly(), getrandpoly()
            f = (1 + a * rand(0:1)) * (1 + b * rand(0:1)) * (1 + c * rand(0:1))
            xs = ParamPunPam.distinct_nonzero_points(ground, max(degree(f) + 1, 1))
            ys = map(x -> evaluate(f, x), xs)
            @test f == Nemo.interpolate(R, xs, ys) == ParamPunPam.fastpolyinterpolate(R, xs, ys)
        end
    end
end

@testset "Transposed Vandermonde solve" begin
    for ground in FIELDS
        R, (x,) = polynomial_ring(ground, ["x"])
        Runiv, xuniv = polynomial_ring(ground, "x")

        # for t == T
        # so that the number of points equals the number of terms
        cases = [
            R(4),
            x,
            x^2,
            2x^3,
            3x^4,
            x + 2,
            x^2 + 1,
            (x - 1) * (x + 4) * (x - 8),
            (x^4 + 4) * (x^2 + 3) * (x + 2),
            (x - 2)^10 * (x^4 + 4)^5 * (x + 11),
            (x - 1)^10 * (x - 2)^20,
            (x - 1)^50 * (x - 2)^70,
            x^20 + 1,
            7(x - 1)^100
        ]
        for f in cases
            ω = [ParamPunPam.randomgenerator(ground)]
            vs = map(m -> evaluate(m, ω), collect(monomials(f)))
            ys = map(i -> evaluate(f, ω .^ i), 0:(length(vs) - 1))
            @test collect(coefficients(f)) == ParamPunPam.solve_transposed_vandermonde(Runiv, vs, ys)
        end
    end
end
