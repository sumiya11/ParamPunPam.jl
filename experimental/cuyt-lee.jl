using Nemo, Random

begin
    using BenchmarkTools, Logging
    import Nemo, Profile
    evalfrac(f, x) = evaluate(numerator(f), x) // evaluate(denominator(f), x)

    macro my_profview(ex)
        esc(:((VSCodeServer.Profile).clear();
        VSCodeServer.Profile.init(n=10^8, delay=0.0001);
        VSCodeServer.Profile.start_timer();
        $ex;
        VSCodeServer.Profile.stop_timer();
        VSCodeServer.view_profile(;)))
    end

    macro my_profview_allocs(ex)
        esc(:((VSCodeServer.Profile).clear();
        VSCodeServer.Profile.Allocs.start(sample_rate=1.0);
        try
            $ex
        finally
            VSCodeServer.Profile.Allocs.stop()
        end;
        VSCodeServer.view_profile_allocs(;)))
    end
end

K = GF(2^62 + 135)
R, x = polynomial_ring(GF(2^62 + 135), ["x$i" for i in 1:15])

case =
    (2x[1] + 3x[2]^2 + 4sum(x) + 1) //
    (5x[1]^3 + 6x[2]^3 + 7sum(x) - 8(sum(x[1:10]) + 3) + 11)
rational_interpolator = ParamPunPam.CuytLee

for _ in 1:30
    n = nvars(R)
    Nt, Dt = 1, 1
    Nd, Dd = total_degree(numerator(case)), total_degree(denominator(case))
    Nd, Dd = Nd < 0 ? 0 : Nd, Dd < 0 ? 0 : Dd
    Nds, Dds = repeat([Nd], n), repeat([Dd], n)
    all_interpolated = false
    interpolator = rational_interpolator(R, Nd, Dd, Nds, Dds, Nt, Dt)
    P, Q = R(0), R(1)
    i = 1
    xs = nothing
    begin
        while !all_interpolated
            all_interpolated = true
            xs = ParamPunPam.get_evaluation_points!(interpolator)
            ys = map(x -> evalfrac(case, x), xs)
            P, Q = ParamPunPam.interpolate!(interpolator, ys)
            dp, dq = total_degree(P), total_degree(Q)
            if dp < total_degree(numerator(case)) || dq < total_degree(denominator(case))
                all_interpolated = false
            end
            if length(P) < length(numerator(case)) || length(Q) < length(denominator(case))
                all_interpolated = false
            end
            @info "" length(xs)
        end
    end
    @info "" length(xs)
    @assert P // Q == case
end
