using Pkg;
Pkg.activate(@__DIR__)

using StructuralIdentifiability
using RationalFunctionFields
using ParamPunPam
using Nemo
using Groebner
using Chairmarks
using PrettyTables
using DataStructures

include("gizmos.jl")

function problems_simplification_identifiable_funcs(;with_states=true)
    # get `benchmarks`
    include(joinpath(dirname(dirname(pathof(StructuralIdentifiability))), "benchmarking", "benchmarks.jl"))
    global benchmarks
    odes = [
        benchmarks[:SIWR_orig],
        benchmarks[:LV],
        benchmarks[:Pharm],
        benchmarks[:SEAIJRC],
        benchmarks[:MAPK_5out],
        benchmarks[:MAPK_6out],
        benchmarks[:Goodwin],
        benchmarks[:HIV],
        benchmarks[:SIRC_forced],
        benchmarks[:Akt],
        benchmarks[:CRN],
        benchmarks[:CD8],
        benchmarks[:QWWC],
        benchmarks[:LLW],
        benchmarks[:HIV2],
        benchmarks[:Treatment],
        benchmarks[:TumorHu],
        benchmarks[:TumorPillis],
        benchmarks[:Biohydrogenation],
        benchmarks[:SLIQR],
        benchmarks[:St],
        benchmarks[:Bilirubin]
    ]
    problems = []
    for ode in odes
        name = string(ode[:name], " with_states=$(with_states)")
        ode = ode[:ode]
        id_funcs =
            StructuralIdentifiability.initial_identifiable_functions(ode, prob_threshold=0.99, with_states=with_states)[1]
        problem = Dict(:name => name, :funcs => id_funcs)
        push!(problems, problem)
    end
    problems
end

function problems_simplification_power_sums()
    problems = []
    for n in 4:7
        for k in (n - 1):(n + 1)
            R, x = polynomial_ring(Nemo.QQ, :x => (1:n,))
            funcs = [power_sums(x, k)]
            name = "power_sums_$(n)_$(k)"
            problem = Dict(:name => name, :funcs => funcs)
            push!(problems, problem)
        end
    end
    problems
end

function a_single_run_of_implementation(funcs)
    rff = RationalFunctionFields.RationalFunctionField(funcs);
    D = (typemax(Int), typemax(Int))
    t = @elapsed RationalFunctionFields.groebner_basis_coeffs(rff, ordering=DegRevLex(), up_to_degree=D)
    @assert length(rff.mqs.cached_groebner_bases) == 1
    gb = first(values(rff.mqs.cached_groebner_bases))
    n_spec, n_reduce = rff.mqs.stats.n_spec_mod_p, rff.mqs.stats.n_red_mod_p

    @info "Time:" t
    @info "Number of evaluations:" n_spec n_reduce
    @info "Interpolated max. degree:" max_degree_of_coeffs(gb)
    @info "Interpolated max. number of terms:" max_terms_of_coeffs(gb)
    return Dict(
        :n_spec=>n_spec,
        :n_red=>n_reduce,
        :our_deg=>max_degree_of_coeffs(gb),
        :terms=>max_terms_of_coeffs(gb),
        :time=>t
    )
end

function compute_max_total_degrees(funcs)
    rff = RationalFunctionFields.RationalFunctionField(funcs);
    degs = paramgb_only_degrees(rff.mqs, ordering=DegRevLex())
    N, D = 0, 0
    for i in degs
        for j in i
            N = max(N, j[1])
            D = max(D, j[2])
        end
    end
    return Dict{Any,Any}(
        :max_deg => (N, D)
    )
end

function try_to_compute_full_gb(funcs)
    rff = RationalFunctionFields.RationalFunctionField(funcs);
    t = @elapsed gb = paramgb(rff.mqs)
    n_spec, n_reduce = rff.mqs.stats.n_spec_mod_p, rff.mqs.stats.n_red_mod_p
    return Dict(
        :n_spec_full=>n_spec,
        :n_red_full=>n_reduce,
        :terms_full=>max_terms_of_coeffs(gb),
        :time_full=>t
    )
end

function compute_sharp_n_spec(funcs; deg=nothing)
    rff = RationalFunctionFields.RationalFunctionField(funcs);
    t = @elapsed gb = paramgb(rff.mqs, up_to_degree=deg)
    n_spec, n_reduce = rff.mqs.stats.n_spec_mod_p, rff.mqs.stats.n_red_mod_p
    return Dict(
        :n_spec_sharp=>n_spec
    )
end

function compute_time_per_eval(funcs)
    p = Nemo.Native.GF(2^62 + 135)
    rff = RationalFunctionFields.RationalFunctionField(funcs);
    StructuralIdentifiability.ParamPunPam.reduce_mod_p!(rff.mqs, p)
    point = rand(p, length(Nemo.gens(StructuralIdentifiability.ParamPunPam.parent_params(rff.mqs))))
    eqs = StructuralIdentifiability.ParamPunPam.specialize_mod_p(rff.mqs, point)
    print("groebner           ")
    b1 = @b Groebner.groebner($eqs, ordering=Groebner.DegRevLex())
    display(b1)
    print("groebner_apply!    ")
    trace, _ = Groebner.groebner_learn(eqs, ordering=DegRevLex())
    b2 = @b Groebner.groebner_apply!($trace, $eqs, ordering=Groebner.DegRevLex())
    display(b2)
    return Dict{Any,Any}(
        :time_per_gb => b1.time,
        :time_per_apply => b2.time
    )
end

function basic_stats(funcs)
    n = length(gens(parent(funcs[1][1])))
    N, D = 0, 0
    for arr in funcs
        D = max(D, total_degree(arr[1]))
        for f in arr[2:end]
            N = max(N, total_degree(f))
        end
    end
    return Dict{Any,Any}(
        :input_n => n,
        :input_deg => (N, D)
    )
end
 
function push_to_table!(table, res)
    for k in keys(table)
        if k in keys(res)
            push!(table[k], res[k])
        else
            push!(table[k], nothing)
        end
    end
end

problems = vcat(
    # problems_simplification_power_sums(),
    problems_simplification_identifiable_funcs(with_states=true),
    problems_simplification_identifiable_funcs(with_states=false),
)

table = OrderedDict(:name => [], :input_n => [], :input_deg => [], :n_spec => [], :n_red => [], :our_deg => [], :max_deg => [], :terms => [], :time => [], :time_per_gb => [], :time_per_apply => [], :terms_full => [], :n_spec_full => [], :n_red_full => [], :time_full => [], :n_spec_sharp => [])

for problem in problems
    name = problem[:name]
    funcs = problem[:funcs]

    println("============= $name ================")
    res0 = basic_stats(funcs)
    res1 = a_single_run_of_implementation(funcs)
    res2 = compute_max_total_degrees(funcs)
    res3 = compute_time_per_eval(funcs)
    # res4 = try_to_compute_full_gb(funcs)
    res4 = Dict()
    res5 = compute_sharp_n_spec(funcs, deg = max.(1, res1[:our_deg]))
    res = merge(res0, res1, res2, res3, res4, res5)
    res[:name] = name
    push_to_table!(table, res)
    println("====================================")
end

pretty_table(reduce(hcat, [v for (k, v) in table]), column_labels=string.(collect(keys(table))))

#=
DEFAULT RUN

┌─────────────────────────────────────────────┬─────────┬───────────┬────────┬───────┬─────────┬──────────┬──────────┬───────────┬─────────────┬────────────────┐
│                                        name │ input_n │ input_deg │ n_spec │ n_red │ our_deg │  max_deg │    terms │      time │ time_per_gb │ time_per_apply │
├─────────────────────────────────────────────┼─────────┼───────────┼────────┼───────┼─────────┼──────────┼──────────┼───────────┼─────────────┼────────────────┤
│              SIWR original with_states=true │      11 │   (10, 0) │     18 │     1 │  (1, 0) │   (1, 0) │   (1, 1) │   2.40801 │   0.0008305 │        0.00043 │
│    Modified LV for testing with_states=true │       6 │    (5, 0) │     30 │     1 │  (2, 0) │   (2, 0) │   (2, 1) │  0.005589 │     6.23e-5 │        3.24e-5 │
│         Goodwin oscillator with_states=true │      11 │   (12, 7) │     62 │     1 │  (2, 2) │   (3, 3) │   (3, 2) │  0.453688 │   0.0004917 │      0.0001526 │
│                        HIV with_states=true │      15 │   (10, 1) │     58 │     1 │  (3, 3) │   (3, 3) │   (1, 1) │  0.692269 │   0.0004114 │        9.87e-5 │
│                SIRS forced with_states=true │      11 │   (17, 0) │     26 │     1 │  (2, 2) │   (2, 2) │   (1, 1) │ 0.0334116 │   0.0011353 │      0.0003902 │
│                Akt pathway with_states=true │      25 │    (7, 1) │     62 │     1 │  (2, 2) │   (2, 2) │   (3, 1) │   5.79549 │   0.0004707 │       0.000152 │
│  Chemical reaction network with_states=true │      12 │    (9, 0) │     18 │     1 │  (1, 0) │   (1, 0) │   (1, 1) │ 0.0222071 │   0.0004292 │      0.0001125 │
│                      SLIQR with_states=true │      10 │   (12, 0) │   1228 │     1 │  (7, 6) │   (7, 6) │ (23, 11) │  0.783162 │   0.0006723 │      0.0001341 │
│                         St with_states=true │      12 │    (7, 2) │    772 │     1 │  (4, 4) │ (11, 11) │ (29, 13) │   11.4842 │   0.0477545 │      0.0094507 │
│              Bilirubin2_io with_states=true │      11 │    (5, 0) │    644 │     1 │  (4, 2) │   (4, 2) │  (18, 6) │  0.904304 │   0.0018434 │      0.0004575 │
│             SIWR original with_states=false │       7 │  (24, 11) │     24 │     1 │  (1, 2) │   (1, 2) │   (1, 1) │   7.49199 │    0.442858 │      0.0696024 │
│   Modified LV for testing with_states=false │       4 │    (2, 3) │     30 │     1 │  (2, 0) │   (2, 0) │   (2, 1) │ 0.0062482 │     4.66e-5 │        2.71e-5 │
│        Goodwin oscillator with_states=false │       7 │   (14, 6) │     34 │     1 │  (2, 1) │   (2, 1) │   (2, 1) │ 0.0325746 │    0.000645 │      0.0001935 │
│                       HIV with_states=false │      10 │    (8, 3) │     56 │     1 │  (3, 2) │   (3, 2) │   (1, 1) │ 0.0604035 │     0.00062 │      0.0001358 │
│               SIRS forced with_states=false │       6 │   (15, 5) │     22 │     1 │  (2, 0) │   (2, 0) │   (1, 1) │   1.48899 │   0.0650385 │      0.0125389 │
│               Akt pathway with_states=false │      16 │  (15, 12) │     46 │     1 │  (1, 1) │   (4, 7) │   (3, 1) │  0.141464 │   0.0018867 │      0.0005381 │
│ Chemical reaction network with_states=false │       6 │    (5, 4) │     34 │     1 │  (1, 2) │   (1, 2) │   (1, 2) │ 0.0314471 │   0.0003445 │      0.0001182 │
│                     SLIQR with_states=false │       6 │    (6, 4) │    716 │     1 │  (7, 6) │   (7, 6) │ (14, 11) │  0.207696 │   0.0003166 │      0.0001008 │
│                        St with_states=false │       9 │  (18, 12) │    772 │     1 │  (4, 4) │  (11, 9) │  (26, 8) │   33.2523 │    0.336453 │      0.0123608 │
│             Bilirubin2_io with_states=false │       7 │    (4, 3) │    316 │     1 │  (4, 2) │   (4, 2) │  (12, 6) │  0.186747 │   0.0003214 │      0.0001259 │
└─────────────────────────────────────────────┴─────────┴───────────┴────────┴───────┴─────────┴──────────┴──────────┴───────────┴─────────────┴────────────────┘

=#
