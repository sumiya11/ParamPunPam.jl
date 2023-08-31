
# Progress bars
const _progressbar_color = :cyan
const _progressbar_value_color = :cyan # :light_grey
progressbar_enabled() =
    Logging.Info <= Logging.min_enabled_level(current_logger()) < Logging.Warn

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
