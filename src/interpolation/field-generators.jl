# Returns a random generator of the multiplicative group of field K
function randomgenerator(K)
    ord = BigInt(order(K) - 1)
    factors = Primes.factor(Vector, ord)
    g = rand(K)
    i = 0
    while !generates_mult_group(ord, factors, g)
        i += 1
        g = rand(K)
        i > ord && error("The characteristic of the base field is too small, sorry")
    end
    g
end

function distinct_generators(k, n)
    ch = Int(characteristic(k))
    ord = ch - 1
    factors = Primes.factor(Vector, ord)
    v = zeros(k, n)
    jj = 2
    for i in 1:n
        while jj < ch
            # kj = k(jj)
            kj = rand(k)
            if generates_mult_group(ord, factors, kj)
                v[i] = kj
                jj += 1
                break
            end
            jj += 1
        end
        jj >= ch && error("The characteristic of the base field is too small, sorry")
    end
    v
end

# Returns alphanp2, a generator of K, such that
# alphanp2/alphanp1 is itself a generator of K,
# and alphanp2 is not in alphan
function distinct_generator(K, alphan, alphanp1)
    ord = Int(order(K) - 1)
    factors = Primes.factor(Vector, ord)
    i = 0
    while true
        i += 1
        i > ord && error("Could not find a distinct generator")
        gamma = rand(K)
        alphanp2 = gamma*alphanp1 
        alphanp2 in alphan && continue
        if generates_mult_group(ord, factors, gamma)
            return alphanp2
        end
    end
    return one(K)
end

function generates_mult_group(ord, factors, kj)
    iszero(kj) && return false
    for p in factors
        d = div(ord, p)
        if isone(kj^d)
            return false
        end
    end
    true
end

# Returns a primitive r-th root of unity in K.
function root_of_unity(K, r)
    ord = order(K)
    factors = Primes.factor(Vector, ord)
    while true
        a = rand(K)
        if generates_mult_group(ord, factors, a)
            return a^div(ord, r)
        end
    end
    # bad branch
    return one(K)
end

# Returns a primitive root of unity in K
# of order >= r and the order itself
function approx_root_of_unity(K, r)
    ord = Int(order(K)) - 1
    @assert r <= ord
    g = first(distinct_generators(K, 1))
    r == ord && return (ord, g)
    factors = Primes.factor(Vector, ord)
    i = 1
    ord1 = ord
    while ord1 >= r && i <= length(factors)
        @assert ord1 >= r
        ord1 = div(ord1, factors[i])
        i += 1
    end
    ord1 = ord1*factors[i - 1]
    @assert ord1 >= r
    (ord1, g^div(ord, ord1))
end
