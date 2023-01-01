
function buchberger(F)
    F = unique(filter(!iszero, F))
    _buchberger(F)
end

function _buchberger(F)
    G = deepcopy(F)
    Spolys = [spoly(f, g) for f in G, g in G if f != g]
    while !isempty(Spolys)
        f = first(Spolys)
        Spolys = Spolys[2:end]
        _, r = divrem(f, G)
        if !iszero(r)
            append!(Spolys, [spoly(f, g) for g in G])
            push!(G, f)
        end
    end
    G
end

function spoly(f, g)
    lcmfg = lcm(leading_term(f), leading_term(g))
    f*divexact(lcmfg, leading_term(f)) - g*divexact(lcmfg, leading_term(g))
end
