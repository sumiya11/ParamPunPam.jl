
const _random_number_bound = 2^10 - 1

random_point(gf::Nemo.fpField) = rand(gf)
random_point(gff::Nemo.FpField) = rand(gff)
random_point(fqn::Nemo.fqPolyRepField) = rand(fqn)
random_point(fqn::Nemo.FqField) = rand(fqn)
random_point(::Nemo.ZZRing) = ZZ(rand(1:_random_number_bound))
random_point(::Nemo.QQField) = QQ(random_point(ZZ))

function random_point(ring)
    K = base_ring(ring)
    map(_ -> random_point(K), 1:nvars(ring))
end

# returns an array of L distinct points from the given field
function distinct_nonzero_points(field, L, prev=nothing)
    @label Start
    ans = [random_point(field) for _ in 1:L]
    if prev !== nothing
        ans = vcat(prev, ans)
    end
    if any(iszero, ans)
        @goto Start
    end
    if length(unique(ans)) != length(ans)
        @goto Start
    end
    ans
end

homogenize(ring; varname="x0") = first(
    polynomial_ring(
        base_ring(ring),
        append!([varname], map(string, AbstractAlgebra.symbols(ring))),
        internal_ordering=Nemo.internal_ordering(ring)
    )
)
dehomogenize(ring) = first(
    polynomial_ring(
        base_ring(ring),
        map(string, AbstractAlgebra.symbols(ring))[2:end],
        internal_ordering=Nemo.internal_ordering(ring)
    )
)

univariatize(::Type{Ring}, ring; varname="x") where {Ring <: AbstractAlgebra.MPolyRing} =
    first(polynomial_ring(base_ring(ring), [varname]))
univariatize(::Type{Ring}, ring; varname="x") where {Ring <: AbstractAlgebra.PolyRing} =
    first(polynomial_ring(base_ring(ring), varname))

function getboundsinfo(f)
    (
        totaldeg=max(0, total_degree(f)),
        nterms=max(1, length(f)),
        partialdegs=map(x -> max(0, degree(f, x)), gens(parent(f)))
    )
end

function getboundsinfo(f::AbstractAlgebra.Generic.FracFieldElem)
    nummi = getboundsinfo(numerator(f))
    denmi = getboundsinfo(denominator(f))
    (
        numtotaldeg=nummi.totaldeg,
        numnterms=nummi.nterms,
        numpartialdegs=nummi.partialdegs,
        dentotaldeg=denmi.totaldeg,
        dennterms=denmi.nterms,
        denpartialdegs=denmi.partialdegs
    )
end

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
