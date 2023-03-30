
const _random_number_bound = 2^10-1

random_point(gf::Nemo.GaloisField) = rand(gf)
random_point(gff::Nemo.GaloisFmpzField) = rand(gff)
random_point(fqn::Nemo.FqNmodFiniteField) = rand(fqn)
random_point(fqn::Nemo.FqFiniteField) = rand(fqn)
random_point(::Nemo.FlintIntegerRing) = ZZ(rand(1:_random_number_bound))
random_point(::Nemo.FlintRationalField) = QQ(random_point(ZZ))

function random_point(ring)
    K = base_ring(ring)
    map(_ -> random_point(K), 1:nvars(ring))
end

# returns an array of L distinct points from the given field
function distinct_points(field, L)
    ans = [random_point(field) for _ in 1:L]
    allunique(ans) && return ans
    distinct_points(field, L)
end

homogenize(ring; varname="x0") = first(PolynomialRing(base_ring(ring), append!([varname], map(string, AbstractAlgebra.symbols(ring)))))
dehomogenize(ring) = first(PolynomialRing(base_ring(ring), map(string, AbstractAlgebra.symbols(ring))[2:end]))

univariatize(::Type{Ring}, ring; varname="x") where {Ring<:AbstractAlgebra.MPolyRing} = first(PolynomialRing(base_ring(ring), [varname]))
univariatize(::Type{Ring}, ring; varname="x") where {Ring<:AbstractAlgebra.PolyRing} = first(PolynomialRing(base_ring(ring), varname))

function getboundsinfo(f)
    (
        totaldeg=max(0, total_degree(f)),
        nterms=max(1, length(f)),
        partialdegs=map(x -> max(0, degree(f, x)), gens(parent(f)))
    )
end

function getboundsinfo(f::AbstractAlgebra.Generic.Frac)
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
