
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
