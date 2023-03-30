
LV_gens = ParamPunPam.load_generators("C:\\data\\projects\\ParamPunPam.jl\\examples\\LV_gens.txt")

I = ParamPunPam.generators_to_saturated_ideal(LV_gens)

ParamPunPam.paramgb(I)

# 

goodwin = ParamPunPam.load_generators("C:\\data\\projects\\ParamPunPam.jl\\examples\\Goodwin.txt")

I = ParamPunPam.generators_to_saturated_ideal(goodwin)

goodwin_gb = ParamPunPam.paramgb(I)

# 

sontag = ParamPunPam.load_generators("C:\\data\\projects\\ParamPunPam.jl\\examples\\Sontag_gen.txt");

I = ParamPunPam.generators_to_saturated_ideal(sontag);

sontag_gb = ParamPunPam.paramgb(I)


