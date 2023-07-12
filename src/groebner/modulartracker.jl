
mutable struct ModularTracker{F}
    # Current finite field
    ff::F
    # The product of all used prime numbers
    primeproduct::BigInt

    function ModularTracker(blackbox)
        ff = Nemo.GF(2^62 + 135)
        # ff = Nemo.GF(2^31 - 1)
        new{typeof(ff)}(ff, BigInt(1))
    end
end

function randluckypoint(modular::ModularTracker)
    rand(modular.ff)
end

