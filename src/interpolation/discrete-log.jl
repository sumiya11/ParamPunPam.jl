# Discrete logarithms in F_q^n

# Stores the field K and some information about it
# to speed up discrete logarithms
mutable struct PrecomputedField{Field}
    K::Field
    ordmult::Int
    factors::Vector{Pair{Int, Int}}
    extensiondeg::Int

    function PrecomputedField(K::Field) where {Field}
        Nemo.order(K) > typemax(Int) && @warn "The field is too large for discrete logarithms."
        ordmult = Int(Nemo.order(K) - 1)
        factors = collect(Primes.factor(Dict, ordmult))
        new{Field}(K, ordmult, factors, Nemo.degree(K))
    end
end

# Stores preallocated buffers to speed discrete logarithms 
mutable struct DiscreteLogBuffers{Field, I}
    PF::PrecomputedField{Field}
    xibuf::Vector{Int}
    pibuf::Vector{Int}
    baby::Dict{I, Int}

    function DiscreteLogBuffers(PF::PrecomputedField{Field}) where {Field}
        nfactors = length(PF.factors)
        largest = maximum(f -> first(f), PF.factors)
        baby = sizehint!(Dict{elem_type(PF.K), Int}(), isqrt(largest))
        new{Field, elem_type(PF.K)}(
            PF, 
            Vector{Int}(undef, nfactors), 
            Vector{Int}(undef, nfactors), 
            baby
        )
    end
end

# Dispatch between direct and baby-giant discrete log algorithms 
# is based on the size of the base order.
# If the order is < _threshold_direct_case(), direct algorithm is used.
# Based on profile results in perf/discrete-logs.ipynb
_threshold_direct_case() = 2^5

# Solves a^x = y (mod p) for x
# ord is the order of a in Z/Zp, factors is an array of factors of p-1
discrete_log(a::I, y::I, buf; ord_isprime=false) where {I} = 
    discrete_log(a, y, buf.PF.ordmult, buf; ord_isprime=ord_isprime)

function discrete_log(a::I, y::I, ord::T, buf; ord_isprime=false) where {I, T}
    if ord < _threshold_direct_case()
        direct_discrete_log(a, y, ord, buf)
    else
        if !ord_isprime
            # make sure pohlig_hellman_discrete_log 
            # can not call itself recursively
            pohlig_hellman_discrete_log(a, y, ord, buf)
        else
            babystep_giantstep_discrete_log(a, y, ord, buf)
        end
    end
end

# Solves a^x = y (mod p) for x using the Pohlig-Hellman algorithm.
# 
# ord is the order of a in Z/Zp.
# factors is a dictionary of prime factors of ord (with multiplicities).
# xibuf and pibuf are buffers used to store intermediate results.
function pohlig_hellman_discrete_log(a::I, y::I, ord::T, buf) where {I, T}
    PF = buf.PF
    @inbounds for i in 1:length(PF.factors)
        (pi, di) = PF.factors[i]
        ai, yi = a, y
        xi = zero(T)
        cc = one(ai)
        for j in 0:di-1
            pij = pi^j
            cij = div(ord, pij*pi)
            aij = ai^(cij*pij)
            yij = (yi*inv(cc))^cij
            xij = discrete_log(aij, yij, pi, buf, ord_isprime=true)
            tij = xij*pij
            cc *= ai^tij
            xi += tij
        end
        buf.xibuf[i] = xi
        buf.pibuf[i] = pi^di
    end
    mod(Nemo.crt(buf.xibuf, buf.pibuf), ord)
end

# Solves a^x = y (mod p) for x using the baby-step giant-step algorithm.
# ord is the order of a in Z/Zp.
function babystep_giantstep_discrete_log(a::I, y::I, ord::T, buf) where {I, T}
    # the size of the giant step
    m = Int(isqrt(ord) + 1)
    # baby steps
    baby = buf.baby
    # this does nothing if the size of baby is enough
    sizehint!(baby, m)
    ai = one(a)
    for i in 0:m-1
        baby[ai] = i
        ai *= a
    end
    # find a match
    ainvm = inv(a)^m
    giantstep = y
    i = 0
    while !haskey(baby, giantstep)
        giantstep *= ainvm
        i += 1
    end
    laststep = baby[giantstep]
    empty!(baby)
    i*m + laststep
end

# Solves a^x = y (mod p).
function direct_discrete_log(a::I, y::I, ord, buf) where {I}
    i = 0
    ai = one(a)
    while ai != y
        ai *= a
        i += 1
    end
    i
end
