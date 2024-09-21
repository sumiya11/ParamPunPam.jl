cases = []

R, (x1, x2) = polynomial_ring(Nemo.Native.GF(2^62 + 135), ["x1", "x2"])
append!(
    cases,
    [
        R(1) // R(1),
        R(1) // R(5),
        R(0) // R(1),
        R(5) // R(9),
        (x1 + x2 + 5) // R(1),
        (x1 + x2 + 5) // R(3),
        x1 // (x1 + 9),
        R(1) // (x1 + 9),
        1 // (x1 * x2),
        1 // (x1 * x2)^3,
        (x1 + 2) // (x1 + 3),
        (x1 + x2)^5 // (x1 + x2 + 5)^8
    ]
)

R, (x0, x1, x2, x3, x4) = polynomial_ring(Nemo.Native.GF(2^62 + 135), ["x0", "x1", "x2", "x3", "x4"])
append!(cases, [(x1 * x4 - x2 * x3) // R(1)])

R, xi = polynomial_ring(Nemo.Native.GF(2^62 + 135), ["x$i" for i in 1:10])
append!(
    cases,
    [
        sum(xi) // (prod(xi)),
        (xi[1] + xi[10]) // (xi[1] - xi[10] - 2),
        (xi[1] + xi[3] + xi[5] + xi[7])^5 // (xi[1] + 2),
        sum(xi) // prod(xi),
        sum(xi) // (prod(xi) + sum(xi) + 2),
        (sum(xi) + 3) // prod(xi),
        (3xi[1] + 2xi[2] + 5xi[3] - 1) // (7xi[4] + 11xi[5] + 13xi[6])
    ]
)

evalfrac(f, x) = evaluate(numerator(f), x) // evaluate(denominator(f), x)

@testset "van-der-Hoeven-Lecerf & Cuyt-Lee" begin
    for rational_interpolator in [ParamPunPam.VanDerHoevenLecerf, ParamPunPam.CuytLee]
        @warn "Testing $rational_interpolator"
        for case in cases
            # Check that direct interpolation works
            @info "" case
            R, n = parent(numerator(case)), nvars(parent(numerator(case)))
            K = base_ring(R)
            Nd, Dd = total_degree(numerator(case)), total_degree(denominator(case))
            Nd, Dd = Nd < 0 ? 0 : Nd, Dd < 0 ? 0 : Dd
            Nds, Dds = repeat([Nd], n), repeat([Dd], n)
            Nt, Dt = length(numerator(case)), length(denominator(case))
            Nt, Dt = iszero(Nt) ? 1 : Nt, iszero(Dt) ? 1 : Dt
            vdhl = rational_interpolator(R, Nd, Dd, Nds, Dds, Nt, Dt)
            xs = ParamPunPam.get_evaluation_points!(vdhl)
            ys = map(x -> evalfrac(case, x), xs)
            success, P, Q = ParamPunPam.interpolate!(vdhl, ys)
            @test success && P // Q == case

            # Check that direct interpolation with small T fails
            if rational_interpolator == ParamPunPam.VanDerHoevenLecerf
                if Nt > 1 && Dt > 1
                    Nt, Dt = 1, 1
                    vdhl = rational_interpolator(R, Nd, Dd, Nds, Dds, Nt, Dt)
                    xs = ParamPunPam.get_evaluation_points!(vdhl)
                    ys = map(x -> evalfrac(case, x), xs)
                    success, P, Q = ParamPunPam.interpolate!(vdhl, ys)
                    @test total_degree(P) < total_degree(numerator(case)) ||
                          total_degree(Q) < total_degree(denominator(case)) ||
                          length(P) < length(numerator(case)) ||
                          length(Q) < length(denominator(case))
                end
            end

            # Check that iterative interpolation works
            for _ in 1:20
                subs = map(K, rand(1:(2^20), n))
                new_x = gens(R) .* subs
                new_case = evaluate(case, new_x)

                Nt, Dt = 1, 1
                Nd, Dd = total_degree(numerator(case)), total_degree(denominator(case))
                Nd, Dd = Nd < 0 ? 0 : Nd, Dd < 0 ? 0 : Dd
                Nds, Dds = repeat([Nd], n), repeat([Dd], n)

                all_interpolated = false
                interpolator = rational_interpolator(R, Nd, Dd, Nds, Dds, Nt, Dt)
                P, Q = R(0), R(1)
                i = 1
                while !all_interpolated
                    all_interpolated = true
                    xs = ParamPunPam.get_evaluation_points!(interpolator)
                    ys = map(x -> evalfrac(new_case, x), xs)
                    success, P, Q = ParamPunPam.interpolate!(interpolator, ys)
                    dp, dq = total_degree(P), total_degree(Q)
                    if dp < total_degree(numerator(new_case)) ||
                       dq < total_degree(denominator(new_case)) ||
                       length(P) < length(numerator(new_case)) ||
                       length(Q) < length(denominator(new_case))
                        all_interpolated = false
                    end
                end
                @test success && P // Q == new_case
            end
        end
    end
end
