
mutable struct ModularTracker{F}
    # Current finite field
    finite_field::F
    # The product of all used prime numbers but the last one
    modulo::BigInt

    used_primes::Vector{UInt64}

    function ModularTracker(blackbox)
        finite_field = Nemo.Native.GF(UInt(big(2)^64 - 59))
        new{typeof(finite_field)}(finite_field, BigInt(1), UInt64[])
    end
end

function find_next_lucky_prime!(mt::ModularTracker)
    current_prime = UInt(BigInt(characteristic(mt.finite_field)))
    next_prime = Primes.prevprime(current_prime - 1)
    push!(mt.used_primes, current_prime)
    mt.finite_field = Nemo.Native.GF(next_prime)
    nothing
end
