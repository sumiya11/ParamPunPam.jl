# Adapted from https://discourse.julialang.org/t/expression-parser/41880/7
# code by Alan R. Rogers, Professor of Anthropology, University of Utah

# Adapted from https://github.com/x3042/Exact-reduction-of-ODE-systems/blob/main/src/myeval.jl
# by A. Demin, HSE

# Adapted from https://github.com/sumiya11/RationalFunctionFields.jl/blob/main/src/myeval.jl
# by A. Demin, HSE

function myeval(e::Union{Expr,Symbol,Number}, map::Dict{Symbol,fmpq_mpoly})
    try
        return _myeval(e, map)
    catch ex
        println("Can't parse \"$e\"")
        rethrow(ex)
    end
end

function _myeval(s::Symbol, map::Dict{Symbol,fmpq_mpoly})
    if haskey(map, s)
        return map[s]
    else
        throw(UndefVarError(s))
    end
end

# Numbers are converted to QQ.
function _myeval(x::Number, map::Dict{Symbol,fmpq_mpoly})
    return QQ(x)
end

# a helper definition for floats
function _myeval(x::Float64, map::Dict{Symbol, fmpq_mpoly})
    result = QQ(0)

    # Getting the result from the string representation in order
    # to avoid approximations caused by the float representation
    s = string(x)
    denom = 1
    extra_num = 1
    if occursin(r"[eE]", s)
        s, exp = split(s, r"[eE]")
        if exp[1] == "-"
            denom = QQ(10)^(-parse(Int, exp))
        else
            extra_num = QQ(10)^(parse(Int, exp))
        end
    end
    frac = split(s, ".")
    if length(frac) == 1
        result = QQ(parse(fmpz, s)) * extra_num // denom
    else 
        result = QQ(parse(fmpz, frac[1] * frac[2])) * extra_num // (denom * 10^(length(frac[2])))
    end
    
    # too verbose for now
    # @warn "a possibility of inexact float conversion" from=x to=result
    return result
end

# To parse an expression, convert the head to a singleton
# type, so that Julia can dispatch on that type.
function _myeval(e::Expr, map::Dict{Symbol,fmpq_mpoly})
    return _myeval(Val(e.head), e.args, map)
end

# Call the function named in args[1]
function _myeval(::Val{:call}, args, map::Dict{Symbol,fmpq_mpoly})
    return _myeval(Val(args[1]), args[2:end], map)
end

# Addition
function _myeval(::Val{:+}, args, map::Dict{Symbol,fmpq_mpoly})
    x = 0
    for arg ∈ args
        x += _myeval(arg, map)
    end
    return x
end

# Subtraction and negation
function _myeval(::Val{:-}, args, map::Dict{Symbol,fmpq_mpoly})
    len = length(args)
    if len == 1
        return -_myeval(args[1], map)
    else
        return _myeval(args[1], map) - _myeval(args[2], map)
    end
end

# Multiplication
function _myeval(::Val{:*}, args, map::Dict{Symbol,fmpq_mpoly})
    x = 1
    for arg ∈ args
        x *= _myeval(arg, map)
    end
    return x
end

# Division
function _myeval(::Val{:/}, args, map::Dict{Symbol,fmpq_mpoly})
    return _myeval(args[1], map) // _myeval(args[2], map)
end

# Exponentiation
function _myeval(::Val{:^}, args, map::Dict{Symbol, fmpq_mpoly})
    if 1 != denominator(_myeval(args[2], map))
        @warn "chto-to strannoe happened"
    end

    return _myeval(args[1], map) ^ Int(numerator(_myeval(args[2], map)))
end

"""
Saturates the given ideal by adding 1 - Q*t to its generators
where t is the new introduced variable

Returns the saturated ideal and the t variable
"""
function saturate(I, Q)
    R = parent(I[1])
    base = base_ring(R)
    strings = map(String, symbols(R))
    parentring, vs = Nemo.PolynomialRing(base, [strings..., "t"], ordering=:degrevlex)

    t = last(vs)
    # @info "" base I parentring
    It = [
        evaluate(f, vs[1:end-1])
        # change_base_ring(base, f, parent=parentring)[1]
        for f in I
    ]
    sat = 1 - Q * t
    push!(It, sat)
    It, t
end

"""
Returns an ideal of polynomials in ys of form
    < Fyi * Qx - Fxi * Qy >
for each Fi in genset
Here Q stands for the lcm of the denominators occuring in genset.
So that generators of the resulting ideal are polynomials,
but not fractions
Saturates the ideal with Q by adding 1 - Qt to its generators
"""
function generators_to_saturated_ideal(genset)

    basepolyring = parent(numerator(first(genset)))
    nvariables = length(gens(basepolyring))
    ground = base_ring(basepolyring)

    dens = map(denominator, genset)
    Q = dens[1]
    for d in dens
        Q = lcm(Q, d)
    end


    Fs = map(numerator ∘ (g -> g * Q), genset)

    # @info "" Q Fs

    ystrings = ["y$i" for i in 1:nvariables]
    yoverx, yoverxvars = Nemo.PolynomialRing(
        basepolyring,
        ystrings,
        ordering=:degrevlex
    )

    Fx = Fs
    Fy = map(F -> change_base_ring(basepolyring, F, parent=yoverx), Fs)
    Qx = Q
    Qy = change_base_ring(basepolyring, Q, parent=yoverx)

    # Gleb: to cancel Fyi and Qy by their gcdi
    # oops
    I = [
        Fyi * Qx - Fxi * Qy
        for (Fyi, Fxi) in zip(Fy, Fx)
    ]

    I, t = saturate(I, Q)

    # (I=I, yoverx=yoverx, basepolyring=basepolyring,
    #     nvariables=nvariables, ground=ground, ystrings=ystrings,
    #     Q=Q, t=t)

    map(f -> map_coefficients(c -> c // basepolyring(1), f), I)
end

"""
    Loads a set of generators from the file by given filepath
    One could check valid file examples in the RFF/data directory
    
    Returns a set of fractions of Nemo polys over Nemo rationals
"""
function load_generators(filepath)
    lines = []
    open(filepath, "r") do inputs
        lines = map(strip, readlines(inputs))
    end

    @info "Loading from $filepath"
    
    strings = map(String, split(lines[1], ", "))

    S, xs = Nemo.PolynomialRing(Nemo.QQ, strings, ordering=:degrevlex)

    mapping = Dict{Symbol, fmpq_mpoly}(
        Symbol(x) => Nemo.gen(S, i)
        for (i, x) in enumerate(strings)
    )
    
    nemoring, = Nemo.PolynomialRing(Nemo.QQ, strings, ordering=:degrevlex)
    
    generators = Nemo.Generic.Frac{fmpq_mpoly}[]

    for line in lines[2:end]
        polystrings = map(Meta.parse ∘ String ∘ strip, split(line, ", "))
        polys = map(f -> myeval(f, mapping), polystrings)
       
        polys = [
            Nemo.isconstant(S(f)) ? nemoring(Nemo.QQ(f)) : change_base_ring(Nemo.QQ, f, parent=nemoring)
            for f in polys
        ]

        append!(generators, [ g // polys[1] for g in polys[2:end] ])
    end
       
    @info "loaded $(length(generators)) generators from $filepath"
    
    return generators
end


