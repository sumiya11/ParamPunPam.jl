
# Keep track of some statistics
const _runtime_data = Dict()

@noinline __throw_unlucky_cancellation() =
    throw(AssertionError("Unlucky cancellation of coefficients in the Groebner basis!"))

@noinline __throw_something_went_wrong(msg) =
    throw(AssertionError("Something went wrong when computing Groebner bases.\n$msg"))

# Progress bars
const _progressbar_color = :cyan
const _progressbar_value_color = :cyan # :light_grey
const _progressbar_spinner = "⌜⌝⌟⌞"
const _is_progressbar_enabled_globally = Ref{Bool}(true)
enable_progressbar(flag::Bool) = _is_progressbar_enabled_globally[] = flag
is_progressbar_enabled() =
    _is_progressbar_enabled_globally[] && Logging.Info <= Logging.min_enabled_level(current_logger()) < Logging.Warn

# Some other utils..
function evaluate_frac(f, x)
    n, d = numerator(f), denominator(f)
    isone(d) && return evaluate(n, x)
    evaluate(n, x) // evaluate(d, x)
end
