using Nemo
using BenchmarkTools

function maybe_add_symbol(subs, val, tag="internal")
    existed = if haskey(subs, val)
        subs[val]
        true
    else
        subs[val] = gensym(tag)
        # push!(subs, (val=val, sym=subexprs[val]))
        subs[val]
        false
    end
    existed, subs[val]
end

function poly_to_slp_1(f)
    ring = parent(f)
    (x,) = gens(ring)
    field = base_ring(ring)
    subs = []
    subexprs = Dict()
    slp = []

    var = x
    _, var_sym = maybe_add_symbol(subexprs, var, "input")
    push!(slp, (input=var,))

    prev_term = 1
    for (i, term) in enumerate(reverse(collect(terms(f))))
        c = coeff(term, 1)
        m = monomial(term, 1)
        d = degree(m, x)
        var = x
        _, c_sym = maybe_add_symbol(subexprs, c, "const")
        _, var_sym = maybe_add_symbol(subexprs, var, "internal")
        if d == 0
            maybe_add_symbol(subexprs, var^d, "const")
        end
        for i in 2:d
            existed, _ = maybe_add_symbol(subexprs, var^i, "internal")
            if !existed
                push!(slp, (arg1=var, arg2=var^(i - 1), op=(*), res=var^i))
            end
        end
        monom_expr = subexprs[var^d]
        _, term_expr = maybe_add_symbol(subexprs, c * var^d, "internal")
        push!(slp, (arg1=c, arg2=var^d, op=(*), res=c * var^d))
        if i == 1
            prev_term = length(slp)
            continue
        end
        @assert i > 1
        poly_tail = slp[prev_term].res + c * var^d
        _, poly_tail_expr = maybe_add_symbol(subexprs, poly_tail, "internal")
        push!(slp, (arg1=slp[prev_term].res, arg2=c * var^d, op=(+), res=poly_tail))
        prev_term = length(slp)
    end
    o = gensym("output")
    subexprs[slp[end].res] = o
    push!(slp, (output=slp[end].res,))
    slp, subexprs
end

function slp_to_expr(slp, subs; prune_constants=true)
    slp_expr = :()
    input = slp[1]
    output = slp[end]
    for i in 2:(length(slp) - 1)
        line = slp[i]
        arg1_sym = subs[line.arg1]
        arg2_sym = subs[line.arg2]
        if prune_constants
            if occursin(string(:const), string(arg1_sym))
                arg1_sym = parse(Int128, string(line.arg1))
            end
            if occursin(string(:const), string(arg2_sym))
                arg2_sym = parse(Int128, string(line.arg2))
            end
        end
        res_sym = subs[line.res]
        op = line.op
        expr = :($res_sym = $op($arg1_sym, $arg2_sym))
        slp_expr = :($slp_expr; $expr)
    end
    slp_expr
end

function slp_to_func(slp, subs)
    input = subs[slp[1].input]
    output = subs[slp[end].output]
    slp_expr = slp_to_expr(slp, subs)
    func_body = quote
        $slp_expr
    end
    func_name = gensym("func")
    func_expr = :(function $func_name($input::T) where {T}
        $func_body
        return $output
    end)
    func_expr
end

R, (x,) = polynomial_ring(QQ, ["x"])
poly = (2x + 1)^80;

slp, subs = poly_to_slp_1(poly);
slp;
@info "" length(slp) length(subs) subs

@btime slp_to_expr(slp, subs);
func_expr = slp_to_func(slp, subs);

func = eval(func_expr);

@code_llvm debuginfo = :none func(100)
@code_warntype func(100)

R_univ, x_univ = QQ["t"]
R, (x,) = polynomial_ring(QQ, ["x"])

poly = (2x + 1)^80;
poly_univ = (2x_univ + 1)^80;

@btime poly_to_slp_1($poly);
@btime slp_to_func($(poly_to_slp_1(poly))...);
s = slp_to_func(poly_to_slp_1(poly)...)
@time func = eval(s)

@time func(1)

func(Int128(1))
evaluate(poly_univ, QQ(1))

@btime func(Int128(1))
@btime evaluate($poly_univ, $(QQ(1)))
