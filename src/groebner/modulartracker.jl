
mutable struct ModularTracker{F}
    # Current finite field
    ff::F
    # The product of all used prime numbers
    modulo::BigInt

    used_primes::Vector{UInt64}

    function ModularTracker(blackbox)
        ff = Nemo.GF(2^62 + 169)
        new{typeof(ff)}(ff, BigInt(1), UInt64[])
    end
end

function find_next_lucky_prime!(mt::ModularTracker)
    current_prime = UInt(BigInt(characteristic(mt.ff)))
    next_prime = Primes.nextprime(current_prime + 1)
    push!(mt.used_primes, current_prime)
    mt.ff = Nemo.GF(next_prime)
    nothing
end
