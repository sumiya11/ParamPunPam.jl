
function rational_reconstruction(a::I, m::I) where {I<:Union{Int, BigInt}}
    a = mod(a, m)
    if a == 0 || m == 0
        return Nemo.QQ(0, 1)
    end
    if m < 0
        m = -m
    end
    if a < 0
        a = m - a
    end
    if a == 1
        return Nemo.QQ(1, 1)
    end
    bnd = sqrt(float(m) / 2)

    U = (I(1), I(0), m)
    V = (I(0), I(1), a)
    while abs(V[3]) >= bnd
        q = div(U[3], V[3])
        T = U .- q .* V
        U = V
        V = T
    end

    t = abs(V[2])
    r = V[3] * sign(V[2])

    # changed from `<= bnd` to `<= m / bnd`
    # we can speed up this !
    # if t <= bnd && gcd(r, t) == 1
    #    return QQ(r, t)
    # end

    return Nemo.QQ(r, t)

    # throw(DomainError(
    #    :($a//$m), "rational reconstruction of $a (mod $m) does not exist"
    # ))

    # return QQ(0, 1)
end
