
mutable struct ModularTracker{T, F}
    # this are needed
    polys::Vector{T}
    # Current finite field
    ff::F
    # The product of all used prime numbers
    primeproduct::BigInt

    function ModularTracker(polys::Vector{T}) where {T}
        # ff = Nemo.GF(2^62 + 135)
        ff = Nemo.GF(2^31 - 1)
        new{T, typeof(ff)}(polys, ff, BigInt(1))
    end
end

function randluckypoint(modular::ModularTracker)
    rand(modular.ff)
end

