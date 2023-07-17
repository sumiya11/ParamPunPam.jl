
cases = []

R, (x1, x2, x3) = PolynomialRing(Nemo.GF(2^62 + 135), ["x1", "x2", "x3"])
append!(
    cases,
    [
        (poly=R(1), primes=true, kron=true),
        (poly=R(8), primes=true, kron=true),
        (poly=2 * x1, primes=true, kron=true),
        (poly=x1^3 + x1, primes=true, kron=true),
        (poly=x1 * x2, primes=true, kron=true),
        (poly=x2 * x3, primes=true, kron=true),
        (poly=x1 * x3, primes=true, kron=true),
        (poly=x1 * x2 * x3, primes=true, kron=true),
        (poly=x1 + 8, primes=true, kron=true),
        (poly=2x1 + 3x3 + 4x2 + 2, primes=true, kron=true),
        (poly=12321(x1 * x2)^5, primes=true, kron=true),
        (poly=x1^9 + 2x2^7 + 3x3^10, primes=true, kron=true)
    ]
)

R, x = PolynomialRing(Nemo.GF(2^62 + 135), ["x$i" for i in 1:10])
append!(
    cases,
    [
        (poly=5sum(x), primes=true, kron=true),
        (poly=prod(x) - 11, primes=true, kron=true),
        (poly=sum(x)^2, primes=true, kron=true),
        (poly=sum(x)^3, primes=true, kron=true),
        # This is the limit for the Primes approach
        (poly=x[1]^10 - x[10]^12 - 2, primes=true, kron=true),
        # This is the limit for the Kronecker approach
        (poly=x[1]^60 - x[10]^60 - 2, primes=false, kron=true),
        (poly=x[1]^9 + 2x[2]^7 + 3x[3]^10, primes=true, kron=true)
    ]
)

R, x = PolynomialRing(Nemo.GF(2^62 + 135), ["x$i" for i in 1:30])
append!(
    cases,
    [
        # This is the limit for the Kronecker approach
        (poly=5sum(x), primes=true, kron=false),
        (poly=sum(x[1:2:end])^2, primes=true, kron=false),
        (poly=(x[1] + 2x[2] + x[3] + 5)^10, primes=true, kron=false),
    ]
)

@testset "Primes Ben-or-Tiwari" begin
    R, (x1, x2, x3) = PolynomialRing(Nemo.GF(2^62 + 135), ["x1", "x2", "x3"])
    poly = 1x1^4 + 2x2*x3*x1^20 + 3x3 + 1
    T, D = 4, 1
    bot = ParamPunPam.PrimesBenOrTiwari(R, T, D)
    ω = ParamPunPam.startingpoint(bot)
    ωs = map(i -> ω .^ i, 0:(2T - 1))
    ys = map(ω -> evaluate(poly, ω), ωs)
    ParamPunPam.interpolate!(bot, ωs, ys)

    for case in cases
        !case.primes && continue

        poly = case.poly
        R = parent(poly)

        T = length(poly)
        D = total_degree(poly)

        bot = ParamPunPam.PrimesBenOrTiwari(R, T, D)
        ω = ParamPunPam.startingpoint(bot)
        ωs = map(i -> ω .^ i, 0:(2T - 1))
        ys = map(ω -> evaluate(poly, ω), ωs)
        @test ParamPunPam.interpolate!(bot, ωs, ys) == poly

        for _ in 1:3
            # Test for the case when D is an upper bound
            T = length(poly)
            D = total_degree(poly)
            D = D + rand(1:5)
            bot = ParamPunPam.PrimesBenOrTiwari(R, T, D)
            ω = ParamPunPam.startingpoint(bot)
            ωs = map(i -> ω .^ i, 0:(2T - 1))
            ys = map(ω -> evaluate(poly, ω), ωs)
            @test ParamPunPam.interpolate!(bot, ωs, ys) == poly

            # Test for the case when T is an upper bound
            T = length(poly)
            D = total_degree(poly)
            T = T * rand(1:5) + rand(1:5)
            bot = ParamPunPam.PrimesBenOrTiwari(R, T, D)
            ω = ParamPunPam.startingpoint(bot)
            ωs = map(i -> ω .^ i, 0:(2T - 1))
            ys = map(ω -> evaluate(poly, ω), ωs)
            @test ParamPunPam.interpolate!(bot, ωs, ys) == poly
        end
    end
end

@testset "Kronecker Ben-or-Tiwari" begin
    for case in cases
        !case.kron && continue

        poly = case.poly
        R = parent(poly)

        T = length(poly)
        D = total_degree(poly)

        bot = ParamPunPam.KronBenOrTiwari(R, T, D)
        ω = ParamPunPam.startingpoint(bot)
        ωs = map(i -> ω .^ i, 0:(2T - 1))
        ys = map(ω -> evaluate(poly, ω), ωs)
        @test ParamPunPam.interpolate!(bot, ωs, ys) == poly

        for _ in 1:3
            # Test for the case when D is an upper bound
            T = length(poly)
            D = total_degree(poly)
            D = D + rand(1:5)
            bot = ParamPunPam.KronBenOrTiwari(R, T, D)
            ω = ParamPunPam.startingpoint(bot)
            ωs = map(i -> ω .^ i, 0:(2T - 1))
            ys = map(ω -> evaluate(poly, ω), ωs)
            @test ParamPunPam.interpolate!(bot, ωs, ys) == poly

            # Test for the case when T is an upper bound
            T = length(poly)
            D = total_degree(poly)
            T = T * rand(1:5) + rand(1:5)
            bot = ParamPunPam.KronBenOrTiwari(R, T, D)
            ω = ParamPunPam.startingpoint(bot)
            ωs = map(i -> ω .^ i, 0:(2T - 1))
            ys = map(ω -> evaluate(poly, ω), ωs)
            @test ParamPunPam.interpolate!(bot, ωs, ys) == poly
        end
    end
end
