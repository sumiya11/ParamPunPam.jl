import Singular
import AbstractAlgebra
using Logging

global_logger(ConsoleLogger(Logging.Warn))
include((@__DIR__) * "/runmodels.jl")

function aa_to_singular(poly)
    #= relic of a kinder past =#
    Rxx = parent(poly)
    Raa = AbstractAlgebra.base_ring(AbstractAlgebra.base_ring(Rxx))
    Rqq = AbstractAlgebra.base_ring(Raa)
    @assert Rqq == Nemo.QQ
    poly_truepoly =
        Nemo.map_coefficients(c -> (@assert isone(denominator(c)); numerator(c)), poly)
    xstrings = map(string, AbstractAlgebra.gens(Rxx))
    ystrings = map(string, AbstractAlgebra.gens(Raa))
    base, _ = Singular.polynomial_ring(Singular.QQ, ystrings, ordering=:lex)
    new_ring, _ =
        Singular.polynomial_ring(base, xstrings, ordering=AbstractAlgebra.ordering(Raa))
    AbstractAlgebra.change_base_ring(
        AbstractAlgebra.base_ring(new_ring),
        poly_truepoly,
        parent=new_ring
    )
end

for thing in ["Simple compartment", "Goodwin", "SIRS forced", "HIV", "SLIQR", "St"]
    system          = ideals[thing]
    singular_system = map(aa_to_singular, system)
    singular_ring   = parent(singular_system[1])
    singular_ideal  = Singular.Ideal(singular_ring, singular_system)

    println("====================\n$thing")
    ##########################
    println("*** Singular.slimgb:")
    singular_slim = @time begin
        Singular.slimgb(singular_ideal, complete_reduction=true)
    end
    # println(AbstractAlgebra.gens(singular_slim))
    ##########################
    try
        println("*** Singular.std:")
        singular_std = @time begin
            Singular.std(singular_ideal, complete_reduction=true)
        end
        # println(AbstractAlgebra.gens(singular_std))
    catch e
        println("Singular failed.")
        println(e)
    end
    ##########################
    println("*** ParamPunPam.paramgb:")
    parampunpam = @time begin
        ParamPunPam.paramgb(system)
    end
    # println(parampunpam)
end
