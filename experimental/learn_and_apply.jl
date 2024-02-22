###########################
# Computing GBs in two phases: Learn and Apply.
# 
# The Learn phase computes a Groebner basis over Z_p and build the graph of
# computations on the way. The graph stores the directions of polynomial
# reductions produced when computing normal forms and S-polynomials.
#
# The Apply phase evaluates the graph on the original input polynomials.

using Nemo, AbstractAlgebra, Primes
using Groebner

function learn_and_apply(polys)
    polys = filter(!iszero, polys)
    @assert !isempty(polys)
    R = parent(first(polys))
    basis_mod_p, graph = nothing, nothing
    FF = Nemo.GF(2^31 - 1)
    @info "Learn phase"
    if !iszero(characteristic(R))
        # Over Z_p, just compute the graph
        if base_ring(R) isa Nemo.FracField
            K = base_ring(base_ring(base_ring(R)))
            @info "Z_p(a...) path"
            point = map(K, Primes.nextprimes(17, nvars(base_ring(base_ring(R)))))
            polys_spec = map(f -> specialize(f, point), polys)
            basis_mod_p, graph = learn!(polys_spec)
        else
            @info "Z_p path"
            basis_mod_p, graph = learn!(polys)
        end
    elseif base_ring(R) in [Nemo.QQ, AbstractAlgebra.QQ]
        # Over QQ, reduce modulo p, and compute the graph
        @info "QQ path"
        @info "Reduction modulo $(characteristic(FF))"
        polys_mod_p = map(f -> reduce_mod_p(f, FF), polys)
        basis_mod_p, graph = learn!(polys_mod_p)
    elseif base_ring(R) isa AbstractAlgebra.Generic.FracFieldElemField
        # Over QQ(a...), specialize, reduce modulo p, and compute the graph
        @info "QQ(a...) path"
        point = Primes.nextprimes(17, nvars(base_ring(base_ring(R))))
        @info "Specialization at $point"
        @info "Reduction modulo $(characteristic(FF))"
        polys_spec = map(f -> specialize(f, point), polys)
        polys_mod_p = map(f -> reduce_mod_p(f, FF), polys_spec)
        basis_mod_p, graph = learn!(polys_mod_p)
    end
    @info "Apply phase"
    @assert isgroebner(basis_mod_p)
    basis, graph = apply!(graph, polys)
    basis, graph
end

reduce_mod_p(f, FF) = map_coefficients(c -> FF(numerator(c)) // FF(denominator(c)), f)
specialize(f, point) = map_coefficients(
    c -> evaluate(numerator(c), point) // evaluate(denominator(c), point),
    f
)

mutable struct ComputationGraph{Poly}
    polys::Vector{Poly}
    # Adjacency list. A node in the graph is a polynomial. A node is represented
    # with an Int, an index into the polys array. 
    adjlist::Dict{Int, Vector{Int}}
    # Edge --> Some info
    metadata::Dict{NTuple{2, Int}, Vector{Any}}
    function ComputationGraph(polys::Vector{Poly}) where {Poly}
        adjlist = Dict{Int, Vector{Int}}()
        metadata = Dict{NTuple{2, Int}, Vector{Any}}()
        new{Poly}(polys, adjlist, metadata)
    end
end

nnodes(graph) = length(graph.polys)
nedges(graph) = sum(length, values(graph.adjlist))
function add_edge!(graph, ij; meta=nothing)
    i, j = ij
    !haskey(graph.adjlist, i) && (graph.adjlist[i] = Int[])
    push!(graph.adjlist[i], j)
    !haskey(graph.metadata, ij) && (graph.metadata[ij] = [])
    push!(graph.metadata[ij], meta)
end
function clean!(graph)
    for (k, v) in graph.adjlist
        graph.adjlist[k] = append!([graph.adjlist[k][1]], unique(graph.adjlist[k][2:end]))
    end
    nothing
end

function spolynomial(f, g)
    mf, mg = leading_monomial(f), leading_monomial(g)
    m = lcm(mf, mg)
    co_mf, co_mg = divexact(m, mf), divexact(m, mg)
    cf, cg = leading_coefficient(f), leading_coefficient(g)
    c = gcd(cf, cg)
    co_cf, co_cg = divexact(cg, c), divexact(cf, c)
    spoly = co_cf * co_mf * f - co_cg * co_mg * g
    spoly
end

function buchberger_criterion(basis, critical_pairs, i, j)
    fi, fj = basis[i], basis[j]
    m = gcd(leading_monomial(fi), leading_monomial(fj))
    isone(m)
end

ext(i, j) = (min(i, j), max(i, j))
function staircase_criterion(basis, critical_pairs, i, j)
    fi, fj = basis[i], basis[j]
    m = lcm(leading_monomial(fi), leading_monomial(fj))
    for k in eachindex(basis)
        k == i || k == j && continue
        !(ext(i, k) in critical_pairs) && !(ext(j, k) in critical_pairs) && continue
        if first(divides(m, leading_monomial(basis[k])))
            return true
        end
    end
    false
end

lcm_deg(i, j, basis) = total_degree(lcm(basis[i], basis[j]))
function select_using_normal_strategy!(basis, critical_pairs)
    _, j = findmin(x -> lcm_deg(x[1], x[2], basis), critical_pairs)
    pair = critical_pairs[j]
    deleteat!(critical_pairs, j)
    pair
end

function update!(critical_pairs, basis, j)
    for i in 1:(j - 1)
        push!(critical_pairs, (i, j))
        # if buchberger_criterion(basis, critical_pairs, i, j)
        #     critical_pairs[end] = (0, 0)
        # end
    end
    # # test existing pairs for redundancy
    # n = length(critical_pairs) - (j - 1)
    # for i in 1:n
    #     k, l = critical_pairs[i]
    #     critical_pairs[n + k] == (0, 0) && continue
    #     critical_pairs[n + l] == (0, 0) && continue
    #     # m = max(lcm_deg(critical_pairs[n + k]..., basis), lcm_deg(critical_pairs[n + l]..., basis))
    #     # if lcm_deg(k, l, basis) > m && first(divides(lcm(basis[k], basis[l]), leading_monomial(basis[j])))
    #     #     critical_pairs[i] = (0, 0)
    #     # end
    # end
    # filter!(x -> !(x == (0, 0)), critical_pairs)
    nothing
end

function normalform(p, V)
    path = Int[]
    q, r = divrem(p, V)
    path = findall(!iszero, q)
    r, path
end

function learn!(polys)
    # The classic Buchberger algorithm, but records the graph of computations
    sort!(polys, by=leading_monomial)
    critical_pairs = Vector{NTuple{2, Int}}()
    basis = empty(polys)
    for poly in polys
        push!(basis, poly)
        update!(critical_pairs, basis, length(basis))
    end
    graph = ComputationGraph(basis)
    d = 0
    while !isempty(critical_pairs)
        (i, j) = select_using_normal_strategy!(basis, critical_pairs)
        buchberger_criterion(basis, critical_pairs, i, j) && continue
        # staircase_criterion(basis, critical_pairs, i, j) && continue
        (iszero(d % 10)) &&
            (@info "Critical pair of degree $(lcm_deg(i, j, basis)), pairs left: $(length(critical_pairs))")
        spoly = spolynomial(basis[i], basis[j])
        nf, path = normalform(spoly, basis)
        iszero(nf) && continue
        push!(basis, nf)
        update!(critical_pairs, basis, length(basis))
        add_edge!(graph, (length(basis), i), meta=(:spoly, i, j))
        for idx in path
            add_edge!(graph, (length(basis), idx), meta=(:reduction,))
        end
        d += 1
    end
    clean!(graph)
    @info "Learned a graph with $(nnodes(graph)) nodes and $(nedges(graph)) edges"
    basis, graph
end

# TODO: reduce lazily
function apply!(graph, polys)
    R = parent(first(polys))
    m = length(polys)
    n = nnodes(graph)
    basis = copy(polys)
    @info "Input: $m polynomials, DParamPunPamred basis: $n polynomials"
    @info "Unwinding the computation graph.."
    metadata = deepcopy(graph.metadata)
    for i in (m + 1):n
        @info "Constructing element $i"
        # Construct the new polynomial from the polynomials basis[1:i-1]
        f = zero(R)
        for j in graph.adjlist[i]
            meta = pop!(metadata[i, j])
            if first(meta) === :spoly
                _, k, l = meta
                f = spolynomial(basis[k], basis[l])
            elseif first(meta) === :reduction
                _ = meta
                g = basis[j]
                _, f = divrem(f, g)
            else
                throw(ArgumentError("Unknown metadata: $meta"))
            end
        end
        push!(basis, f)
    end
    @label Letsgo
    filter!(!iszero, basis)
    sort!(basis, by=leading_monomial)
    updated = false
    for i in 2:length(basis)
        iszero(basis[i]) && continue
        for j in 1:(i - 1)
            iszero(basis[j]) && continue
            q, r = divrem(basis[i], basis[j])
            basis[i] = r
            updated = updated || !iszero(q)
        end
    end
    updated && @goto Letsgo
    basis, graph
end

###########################

# Some tests!
begin
    R, (x, y, z) = Nemo.QQ["x", "y", "z"]
    Ry, (a, b, c) = Nemo.fraction_field(R)["a", "b", "c"]
    @assert repr(reduce_mod_p(x + 7y, Nemo.GF(5))) === "x + 2*y"
    @assert repr(specialize((x) * a + (z // y) * c, Nemo.QQ.([2, 3, 5]))) === "2*a + 5//3*c"

    c = Groebner.cyclicn(4, ground=Nemo.GF(2^31 - 1), ordering=:degrevlex)
    gb, graph = learn!(c)
    @assert isgroebner(gb)
    c = Groebner.katsuran(3, ground=AbstractAlgebra.QQ, ordering=:degrevlex)
    gb, graph = learn!(c)
    @assert isgroebner(gb)
    c = Groebner.cyclicn(4, ground=Nemo.QQ, ordering=:degrevlex)
    gb, graph = learn!(c)
    @assert isgroebner(gb)
end

# # Over QQ
# begin
#     F = Groebner.katsuran(6, ordering=:degrevlex, ground=Nemo.QQ)
#     gb, graph = learn_and_apply(F)
#     # The basis is correct! But the check fails when the basis is not autoreduced
#     # @assert Groebner.isgroebner(gb)
#     gb, graph
# end;

# begin
#     F = Groebner.katsuran(6, ordering=:degrevlex, ground=Nemo.GF(2^31-1))
#     # @benchmark groebner($F)
# end

# Over QQ(a...)
begin
    Rparam, (a, b) = polynomial_ring(Nemo.QQ, ["a", "b"])
    R, (x, y, z) =
        polynomial_ring(Nemo.fraction_field(Rparam), ["x", "y", "z"], ordering=:degrevlex)
    F = [x^2 + x + (a + 1)^5, x * y + b * y * z + 1 // (a * b)^12, x * z + z + b]
    gb, graph = learn_and_apply(F)

    Rparam, (a, b) = polynomial_ring(Nemo.GF(2^31 - 1), ["a", "b"])
    R, (x, y, z) =
        polynomial_ring(Nemo.fraction_field(Rparam), ["x", "y", "z"], ordering=:degrevlex)
    F = [x^2 + x + (a + 1)^5, x * y + b * y * z + 1 // (a * b)^12, x * z + z + b]
    gb, graph = learn_and_apply(F)
    # also correct
end;

# Interpolation approach fails here, as the degrees of parameters approach 15-20.
# begin
#     Rparam, a = polynomial_ring(Nemo.QQ, ["a$i" for i in 1:5])
#     R, x = polynomial_ring(Nemo.fraction_field(Rparam), ["x$i" for i in 1:5], ordering=:degrevlex)
#     F = [
#         x[1]^2 + x[2] + (a[1] + a[2]),
#         (a[3] * a[4] + 1) * x[1] * x[2] + a[2] * x[2] * x[3] + 1 // (a[1] * a[2])^2,
#         (a[1] * a[2] * a[3]) * x[1] * x[3] + x[3] + a[2],
#         x[3] * x[4]^2 + x[2] * x[4] + a[5],
#         x[4] * x[5]^2 + x[5] + a[1],
#     ]
#     # gb, graph = learn_and_apply(F)
# end;

begin
    Ra, (beta, gamma, v, psi) =
        polynomial_ring(Nemo.GF(2^31 - 1), ["beta", "gamma", "v", "psi"])
    R, (t, y1, y2, y3, y4) = polynomial_ring(
        Nemo.fraction_field(Ra),
        ["t", "y1", "y2", "y3", "y4"],
        ordering=:degrevlex
    )

    SEIR_1_io = [
        (-beta * gamma * psi + 2 * beta * gamma + beta * v + beta * psi) * y2^3 * y4^2 +
        (beta * gamma * psi - 2 * beta * gamma - beta * v - beta * psi) * y2^3 * y4 +
        (beta * gamma * psi - 2 * beta * gamma - beta * v - beta * psi) * y2^2 * y3 * y4 +
        (beta * gamma * psi - 2 * beta * gamma - beta * v - beta * psi) * y2^2 * y4^2 +
        (gamma^3 * psi^2 - gamma^3 * psi - gamma^2 * v * psi - gamma^2 * psi^2) *
        y1 *
        y2 *
        y4 +
        (
            -2 * gamma^3 * psi^2 +
            2 * gamma^3 * psi +
            2 * gamma^2 * v * psi +
            2 * gamma^2 * psi^2
        ) *
        y1 *
        y2 +
        (-gamma^3 * psi^2 + gamma^3 * psi + gamma^2 * v * psi + gamma^2 * psi^2) * y1 * y3 +
        (-gamma^3 * psi^2 + gamma^3 * psi + gamma^2 * v * psi + gamma^2 * psi^2) * y1 * y4,
        (-2 * beta * gamma * psi + 4 * beta * gamma + 2 * beta * v + 2 * beta * psi) *
        y2^3 *
        y4^2 +
        (
            2 * gamma^2 * psi^2 - 2 * gamma^2 * psi - 2 * gamma * v * psi -
            2 * gamma * psi^2
        ) *
        y1 *
        y2^2 *
        y4 +
        (2 * beta * gamma * psi - 4 * beta * gamma - 2 * beta * v - 2 * beta * psi) *
        y2^3 *
        y4 +
        (2 * beta * gamma * psi - 4 * beta * gamma - 2 * beta * v - 2 * beta * psi) *
        y2^2 *
        y3 *
        y4 +
        (2 * beta * gamma * psi - 4 * beta * gamma - 2 * beta * v - 2 * beta * psi) *
        y2^2 *
        y4^2 +
        (
            -4 * gamma^2 * psi^2 +
            4 * gamma^2 * psi +
            4 * gamma * v * psi +
            4 * gamma * psi^2
        ) *
        y1 *
        y2^2 +
        (
            -2 * gamma^2 * psi^2 +
            2 * gamma^2 * psi +
            2 * gamma * v * psi +
            2 * gamma * psi^2
        ) *
        y1 *
        y2 *
        y3 +
        (
            -2 * gamma^2 * psi^2 +
            2 * gamma^2 * psi +
            2 * gamma * v * psi +
            2 * gamma * psi^2
        ) *
        y1 *
        y2 *
        y4,
        0,
        y2^3 * y4^2 - y2^3 * y4 - y2^2 * y3 * y4 - y2^2 * y4^2 +
        (-gamma^3 * psi + gamma^3 + gamma^2 * v + gamma^2 * psi) * y4,
        -y2^3 * y4^2 +
        y2^3 * y4 +
        y2^2 * y3 * y4 +
        y2^2 * y4^2 +
        (gamma^3 * psi - gamma^3 - gamma^2 * v - gamma^2 * psi) * y4,
        beta * y2^3 * y4^2 - beta * y2^3 * y4 - beta * y2^2 * y3 * y4 - beta * y2^2 * y4^2 +
        (-gamma * psi^2 + gamma * psi + v * psi + psi^2) * y1 * y2^2,
        y2^3 * y4^2 - y2^3 * y4 - y2^2 * y3 * y4 - y2^2 * y4^2 +
        (-gamma^2 * psi + gamma^2 + gamma * v + gamma * psi) * y2 * y4,
        -y2^3 * y4^2 + y2^3 * y4 + y2^2 * y3 * y4 + y2^2 * y4^2 + gamma^2 * y2 * y4^2 -
        gamma^2 * y2 * y4 - gamma^2 * y3 * y4 - gamma^2 * y4^2,
        (gamma * psi^2 - gamma * psi - v * psi - psi^2) * y1 * y2^3 * y4 +
        (-beta * gamma * psi + 2 * beta * gamma + beta * v + beta * psi) * y2^3 * y4^2 +
        (-2 * gamma * psi^2 + 2 * gamma * psi + 2 * v * psi + 2 * psi^2) * y1 * y2^3 +
        (-gamma * psi^2 + gamma * psi + v * psi + psi^2) * y1 * y2^2 * y3 +
        (-gamma * psi^2 + gamma * psi + v * psi + psi^2) * y1 * y2^2 * y4 +
        (beta * gamma * psi - 2 * beta * gamma - beta * v - beta * psi) * y2^3 * y4 +
        (beta * gamma * psi - 2 * beta * gamma - beta * v - beta * psi) * y2^2 * y3 * y4 +
        (beta * gamma * psi - 2 * beta * gamma - beta * v - beta * psi) * y2^2 * y4^2,
        (
            2 * gamma^2 * psi^2 - 2 * gamma^2 * psi - 2 * gamma * v * psi -
            2 * gamma * psi^2
        ) *
        y1 *
        y2^3 *
        y4 +
        (
            3 * gamma^2 * psi^2 - 3 * gamma^2 * psi - 3 * gamma * v * psi -
            3 * gamma * psi^2
        ) *
        y1 *
        y2^2 *
        y3 *
        y4 +
        (
            -2 * beta * gamma^2 * psi + 2 * beta * gamma^2 - 3 * beta * gamma * v * psi +
            5 * beta * gamma * v +
            2 * beta * gamma * psi +
            3 * beta * v * psi
        ) *
        y2^3 *
        y4^2 +
        (
            -2 * gamma^2 * psi^2 +
            2 * gamma^2 * psi +
            2 * gamma * v * psi +
            2 * gamma * psi^2
        ) *
        y1 *
        y2^3 +
        (
            -5 * gamma^2 * psi^2 +
            5 * gamma^2 * psi +
            5 * gamma * v * psi +
            5 * gamma * psi^2
        ) *
        y1 *
        y2^2 *
        y3 +
        (
            -2 * gamma^2 * psi^2 +
            2 * gamma^2 * psi +
            2 * gamma * v * psi +
            2 * gamma * psi^2
        ) *
        y1 *
        y2^2 *
        y4 +
        (
            2 * beta * gamma^2 * psi - 2 * beta * gamma^2 + 3 * beta * gamma * v * psi -
            5 * beta * gamma * v - 2 * beta * gamma * psi - 3 * beta * v * psi
        ) *
        y2^3 *
        y4 +
        (
            -3 * gamma^2 * psi^2 +
            3 * gamma^2 * psi +
            3 * gamma * v * psi +
            3 * gamma * psi^2
        ) *
        y1 *
        y2 *
        y3 *
        y4 +
        (
            2 * beta * gamma^2 * psi - 2 * beta * gamma^2 + 3 * beta * gamma * v * psi -
            5 * beta * gamma * v - 2 * beta * gamma * psi - 3 * beta * v * psi
        ) *
        y2^2 *
        y3 *
        y4 +
        (
            2 * beta * gamma^2 * psi - 2 * beta * gamma^2 + 3 * beta * gamma * v * psi -
            5 * beta * gamma * v - 2 * beta * gamma * psi - 3 * beta * v * psi
        ) *
        y2^2 *
        y4^2,
        (-gamma * psi + 2 * gamma + v + psi) * y2^3 * y4^2 +
        (gamma * psi - 2 * gamma - v - psi) * y2^3 * y4 +
        (gamma * psi - 2 * gamma - v - psi) * y2^2 * y3 * y4 +
        (gamma^2 * psi - gamma^2 - gamma * v - 2 * gamma - v - psi) * y2^2 * y4^2 +
        (-2 * gamma^2 * psi + 2 * gamma^2 + 2 * gamma * v + 2 * gamma * psi) * y2^2 * y4 +
        (-gamma^2 * psi + gamma^2 + gamma * v + gamma * psi) * y2 * y3 * y4 +
        (-gamma^2 * psi + gamma^2 + gamma * v + gamma * psi) * y2 * y4^2,
        (gamma * psi^2 - gamma * psi - v * psi - psi^2) * y1 * y2^4 * y3 * y4 +
        (-gamma * psi^2 + gamma * psi + v * psi + psi^2) * y1 * y2^4 * y3 +
        (-gamma * psi^2 + gamma * psi + v * psi + psi^2) * y1 * y2^3 * y3 * y4 +
        (-beta * gamma^2 * v * psi + beta * gamma^2 * v + beta * gamma * v * psi) *
        y2^3 *
        y4^2 +
        (beta * gamma^2 * v * psi - beta * gamma^2 * v - beta * gamma * v * psi) *
        y2^3 *
        y4 +
        (beta * gamma^2 * v * psi - beta * gamma^2 * v - beta * gamma * v * psi) *
        y2^2 *
        y3 *
        y4 +
        (beta * gamma^2 * v * psi - beta * gamma^2 * v - beta * gamma * v * psi) *
        y2^2 *
        y4^2,
        (gamma * psi - 2 * gamma - v - psi) * y2^3 * y4^2 +
        (-gamma * psi + 2 * gamma + v + psi) * y2^3 * y4 +
        (-gamma * psi + 2 * gamma + v + psi) * y2^2 * y3 * y4 +
        (-gamma^2 * psi + gamma^2 + gamma * v + 2 * gamma + v + psi) * y2^2 * y4^2 +
        (2 * gamma^2 * psi - 2 * gamma^2 - 2 * gamma * v - 2 * gamma * psi) * y2^2 * y4 +
        (gamma^2 * psi - gamma^2 - gamma * v - gamma * psi) * y2 * y3 * y4 +
        (gamma^2 * psi - gamma^2 - gamma * v - gamma * psi) * y2 * y4^2,
        (gamma * psi - 2 * gamma - v - psi) * y2^3 * y4^2 +
        (-gamma * psi + 2 * gamma + v + psi) * y2^3 * y4 +
        (-gamma * psi + 2 * gamma + v + psi) * y2^2 * y3 * y4 +
        (-gamma * psi + 2 * gamma + v + psi) * y2^2 * y4^2 +
        (-gamma^3 * psi + gamma^3 + gamma^2 * v + gamma^2 * psi) * y2 * y4^2 +
        (2 * gamma^3 * psi - 2 * gamma^3 - 2 * gamma^2 * v - 2 * gamma^2 * psi) * y2 * y4 +
        (gamma^3 * psi - gamma^3 - gamma^2 * v - gamma^2 * psi) * y3 * y4 +
        (gamma^3 * psi - gamma^3 - gamma^2 * v - gamma^2 * psi) * y4^2,
        (
            -beta * gamma^2 * psi + beta * gamma^2 - beta * gamma * v * psi +
            2 * beta * gamma * v +
            beta * gamma * psi +
            beta * v * psi
        ) *
        y2^3 *
        y4^2 +
        (gamma^3 * psi^2 - gamma^3 * psi - gamma^2 * v * psi - gamma^2 * psi^2) *
        y1 *
        y2^2 *
        y4 +
        (
            beta * gamma^2 * psi - beta * gamma^2 + beta * gamma * v * psi -
            2 * beta * gamma * v - beta * gamma * psi - beta * v * psi
        ) *
        y2^3 *
        y4 +
        (gamma^3 * psi^2 - gamma^3 * psi - gamma^2 * v * psi - gamma^2 * psi^2) *
        y1 *
        y2 *
        y3 *
        y4 +
        (
            beta * gamma^2 * psi - beta * gamma^2 + beta * gamma * v * psi -
            2 * beta * gamma * v - beta * gamma * psi - beta * v * psi
        ) *
        y2^2 *
        y3 *
        y4 +
        (
            beta * gamma^2 * psi - beta * gamma^2 + beta * gamma * v * psi -
            2 * beta * gamma * v - beta * gamma * psi - beta * v * psi
        ) *
        y2^2 *
        y4^2 +
        (-gamma^3 * psi^2 + gamma^3 * psi + gamma^2 * v * psi + gamma^2 * psi^2) *
        y1 *
        y2^2 +
        (
            -2 * gamma^3 * psi^2 +
            2 * gamma^3 * psi +
            2 * gamma^2 * v * psi +
            2 * gamma^2 * psi^2
        ) *
        y1 *
        y2 *
        y3 +
        (-gamma^3 * psi^2 + gamma^3 * psi + gamma^2 * v * psi + gamma^2 * psi^2) *
        y1 *
        y2 *
        y4 +
        (-gamma^3 * psi^2 + gamma^3 * psi + gamma^2 * v * psi + gamma^2 * psi^2) *
        y1 *
        y3 *
        y4,
        beta * y2^3 * y4^2 - beta * y2^3 * y4 - beta * y2^2 * y3 * y4 - beta * y2^2 * y4^2 +
        (-gamma^3 * psi^2 + gamma^3 * psi + gamma^2 * v * psi + gamma^2 * psi^2) * y1,
        (gamma * psi^2 - gamma * psi - v * psi - psi^2) * y1 * y2^4 * y4 +
        (3 * gamma * psi^2 - 3 * gamma * psi - 3 * v * psi - 3 * psi^2) *
        y1 *
        y2^3 *
        y3 *
        y4 +
        (-gamma * psi^2 + gamma * psi + v * psi + psi^2) * y1 * y2^4 +
        (-4 * gamma * psi^2 + 4 * gamma * psi + 4 * v * psi + 4 * psi^2) * y1 * y2^3 * y3 +
        (-gamma * psi^2 + gamma * psi + v * psi + psi^2) * y1 * y2^3 * y4 +
        (-3 * gamma * psi^2 + 3 * gamma * psi + 3 * v * psi + 3 * psi^2) *
        y1 *
        y2^2 *
        y3 *
        y4 +
        (
            -beta * gamma^2 * psi + beta * gamma^2 - 3 * beta * gamma * v * psi +
            4 * beta * gamma * v +
            beta * gamma * psi +
            3 * beta * v * psi
        ) *
        y2^3 *
        y4^2 +
        (
            beta * gamma^2 * psi - beta * gamma^2 + 3 * beta * gamma * v * psi -
            4 * beta * gamma * v - beta * gamma * psi - 3 * beta * v * psi
        ) *
        y2^3 *
        y4 +
        (
            beta * gamma^2 * psi - beta * gamma^2 + 3 * beta * gamma * v * psi -
            4 * beta * gamma * v - beta * gamma * psi - 3 * beta * v * psi
        ) *
        y2^2 *
        y3 *
        y4 +
        (
            beta * gamma^2 * psi - beta * gamma^2 + 3 * beta * gamma * v * psi -
            4 * beta * gamma * v - beta * gamma * psi - 3 * beta * v * psi
        ) *
        y2^2 *
        y4^2,
        2 * beta * y2^3 * y4^2 - 2 * beta * y2^3 * y4 - 2 * beta * y2^2 * y3 * y4 -
        2 * beta * y2^2 * y4^2 +
        (
            -2 * gamma^2 * psi^2 +
            2 * gamma^2 * psi +
            2 * gamma * v * psi +
            2 * gamma * psi^2
        ) *
        y1 *
        y2,
        y2^4 * y4^2 - y2^4 * y4 - y2^3 * y3 * y4 +
        (-gamma - 1) * y2^3 * y4^2 +
        gamma * y2^3 * y4 +
        gamma * y2^2 * y3 * y4 +
        gamma * y2^2 * y4^2,
        -t * y2^3 * y4^2 + t * y2^3 * y4 + t * y2^2 * y3 * y4 + t * y2^2 * y4^2 + 1
    ]

    gb, graph = learn_and_apply(SEIR_1_io)
end

# Too large..
# using StructuralIdentifiability
# begin
#     goodwin = @ODEmodel(
#         x1'(t) = -b * x1(t) + 1 / (c + x4(t)),
#         x2'(t) = alpha * x1(t) - beta * x2(t),
#         x3'(t) = gama * x2(t) - delta * x3(t),
#         x4'(t) = sigma * x4(t) * (gama * x2(t) - delta * x3(t)) / x3(t),
#         y(t) = x1(t)
#     )

#     gens = StructuralIdentifiability.ideal_generators(goodwin)
#     # gb, graph = learn_and_apply(gens);

#     # @time apply!(graph, gens)
# end;

# # Some plotting
# begin
#     return 0  # comment this line if you wish to plot!
#     using Plots
#     using GraphRecipes
#     @show graph.adjlist
#     default(size=(1000, 1000))
#     n = maximum(map(maximum, collect(values(graph.adjlist))))
#     adjlist = [get(graph.adjlist, k, Int[]) for k in 1:n]
#     graphplot(
#         adjlist,
#         names="poly " .* string.(1:n),
#         curves=false,
#         nodeshape=:rect,
#         # self_edge_size=0.25, 
#         # method=:tree
#     )
# end
