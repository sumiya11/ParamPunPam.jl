
# Progress bars
const _progressbar_color = :cyan
const _progressbar_value_color = :cyan # :light_grey
progressbar_enabled() =
    Logging.Info <= Logging.min_enabled_level(current_logger()) < Logging.Warn

const _runtime_data = Dict()

function evaluate_frac(f, x)
    n, d = numerator(f), denominator(f)
    isone(d) && return evaluate(n, x)
    evaluate(n, x) // evaluate(d, x)
end

function groebner_loglevel()
    if current_logger().min_level < Logging.Info
        -3
    else
        0
    end
end

function homogenize_many(fs)
    ring = parent(fs[1])
    newring, hom_vars = PolynomialRing(
        base_ring(ring),
        vcat("X0", map(string, gens(ring))),
        ordering=ordering(ring)
    )
    Fs = empty(fs)
    for f in fs
        D = total_degree(f)
        new_f = zero(newring)
        for term in terms(f)
            cf = coeff(term, 1)
            ev = monomial(term, 1)
            d = total_degree(ev)
            new_f += cf * evaluate(ev, hom_vars[2:end]) * hom_vars[1]^(D - d)
        end
        push!(Fs, new_f)
    end
    return Fs
end

function dehomogenize_many(Fs)
    ring = parent(Fs[1])
    newring, dehom_vars = PolynomialRing(
        base_ring(ring),
        map(string, gens(ring)[2:end]),
        ordering=ordering(ring)
    )
    fs = empty(Fs)
    for F in Fs
        f = evaluate(F, vcat(one(newring), dehom_vars))
        push!(fs, f)
    end
    return fs
end
