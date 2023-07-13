using Nemo

function generate_determinant_circulant(K, n::Int)
    R, xs = PolynomialRing(K, ["x$i" for i in 1:n])
    S = MatrixSpace(R, n, n)
    M = zero(S)
    xnto1(i) = i < n ? xs[i] : one(R)
    for i in 1:2n-1
        if i <= n
            for j in 1:n - i + 1
                M[j, j + i - 1] = xnto1(i)
            end
        else
            for j in (i - n + 1):n
                M[j, j - (i - n)] = xnto1(2n - i + 1)
            end
        end
    end
    
    M, det(M)
end

d11 = deepcopy(d);
d12 = deepcopy(d);

K = GF(2^31-1)
@time M, d = generate_determinant_circulant(K, 12);
length(d)

@time fact = map(first, collect(factor(d)));
map(length, fact)

R = parent(d)
xs = gens(R)
sum([length(coeff(d, [1], [i])) for i in 0:12])

