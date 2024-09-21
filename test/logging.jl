using Logging

function version_agnostic_logger(T, level, stream=stdout)
    @static if VERSION >= v"1.7.0"
        T(level)
    else
        T(stream, level)
    end
end

struct MyLogger <: Logging.AbstractLogger
    _logger::Logging.ConsoleLogger
    MyLogger(loglevel) = new(version_agnostic_logger(Logging.ConsoleLogger, loglevel))
end

Logging.min_enabled_level(logger::MyLogger) = Logging.min_enabled_level(logger._logger)
Logging.shouldlog(logger::MyLogger, level, _module, group, id) =
    Logging.shouldlog(logger._logger, level, _module, group, id)
Logging.handle_message(logger, level, message, _module, group, id, file, line; kwargs...) =
    Logging.handle_message(logger._logger, level, message, _module, group, id, file, line; kwargs...)

loggers_to_test = [
    version_agnostic_logger(ConsoleLogger, Logging.Warn),
    version_agnostic_logger(ConsoleLogger, Logging.Info),
    version_agnostic_logger(ConsoleLogger, Logging.Debug),
    version_agnostic_logger(SimpleLogger, Logging.Info),
    version_agnostic_logger(SimpleLogger, Logging.Warn),
    MyLogger(Logging.Debug),
    MyLogger(Logging.Warn)
]

@testset "Generic logging" begin
    Ra, (a, b, c) = polynomial_ring(Nemo.QQ, ["a", "b", "c"])
    Rx, (x, y, z) = polynomial_ring(Nemo.fraction_field(Ra), ["x", "y", "z"], internal_ordering=:degrevlex)
    F = [y + z + 2^10, x + a^2 // (b + c)]

    for logger in loggers_to_test
        global_logger(logger)
        @test ParamPunPam.paramgb(F) == F
    end

    global_logger(MyLogger(Logging.Debug))
    @test_logs ParamPunPam.paramgb(F)

    newlogger = version_agnostic_logger(Logging.ConsoleLogger, Logging.Info)
    global_logger(newlogger)
end
