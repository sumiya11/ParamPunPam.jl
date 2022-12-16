
LV_gens = ParamPanPam.load_generators("C:\\data\\projects\\ParamPanPam.jl\\examples\\LV_gens.txt")

I = ParamPanPam.generators_to_saturated_ideal(LV_gens)

ParamPanPam.paramgb(I)

# 

goodwin = ParamPanPam.load_generators("C:\\data\\projects\\ParamPanPam.jl\\examples\\Goodwin.txt")

I = ParamPanPam.generators_to_saturated_ideal(goodwin)

goodwin_gb = ParamPanPam.paramgb(I)

# 

sontag = ParamPanPam.load_generators("C:\\data\\projects\\ParamPanPam.jl\\examples\\Sontag_gen.txt");

I = ParamPanPam.generators_to_saturated_ideal(sontag);

sontag_gb = ParamPanPam.paramgb(I)


