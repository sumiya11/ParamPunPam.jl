
LV_gens = ParamPunPam.load_generators((@__DIR__)*"/LV_gens.txt")

I = ParamPunPam.generators_to_saturated_ideal(LV_gens)

ParamPunPam.paramgb(I)

# 

goodwin = ParamPunPam.load_generators((@__DIR__)*"/Goodwin.txt")

I = ParamPunPam.generators_to_saturated_ideal(goodwin)

goodwin_gb = ParamPunPam.paramgb(I)

# 

sontag = ParamPunPam.load_generators((@__DIR__)*"/Sontag_gen.txt");

I = ParamPunPam.generators_to_saturated_ideal(sontag);

sontag_gb = ParamPunPam.paramgb(I)
