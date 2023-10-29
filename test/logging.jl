using Logging

struct MyLogger <: Logging.AbstractLogger
    _logger::Logging.ConsoleLogger
    MyLogger(loglevel) = new(Logging.ConsoleLogger(loglevel))
end

Logging.min_enabled_level(logger::MyLogger) = Logging.min_enabled_level(logger._logger)
Logging.shouldlog(logger::MyLogger, level, _module, group, id) =
    Logging.shouldlog(logger._logger, level, _module, group, id)
Logging.handle_message(logger, level, message, _module, group, id, file, line; kwargs...) =
    Logging.handle_message(
        logger._logger,
        level,
        message,
        _module,
        group,
        id,
        file,
        line;
        kwargs...
    )

loggers_to_test = [
    ConsoleLogger(Logging.Warn),
    ConsoleLogger(Logging.Info),
    ConsoleLogger(Logging.Debug),
    SimpleLogger(Logging.Info),
    SimpleLogger(Logging.Warn),
    MyLogger(Logging.Debug),
    MyLogger(Logging.Warn)
]

@testset "Generic logging" begin
    Ra, (a, b, c) = PolynomialRing(Nemo.QQ, ["a", "b", "c"])
    Rx, (x, y, z) =
        PolynomialRing(Nemo.FractionField(Ra), ["x", "y", "z"], ordering=:degrevlex)
    F = [y + z + 2^10, x + a^2 // (b + c)]

    for logger in loggers_to_test
        global_logger(logger)
        @test ParamPunPam.paramgb(F) == F
    end

    global_logger(MyLogger(Logging.Debug))
    @test_logs ParamPunPam.paramgb(F)

    global_logger(Logging.ConsoleLogger(Logging.Info))
end
