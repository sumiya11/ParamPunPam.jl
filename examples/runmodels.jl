using Logging

using AbstractAlgebra
using Nemo

using StructuralIdentifiability
using StructuralIdentifiability: extract_identifiable_functions_raw

using ParamPunPam

logger = SimpleLogger(stdout, Logging.Debug)

models = [
    # Current result: -a12 - a01 - a21, -a12 - a01, a12*a01
    # Simpler version: a21, a12 + a01, a12*a01
    Dict(
        :name => "Simple compartment",
        :ode => @ODEmodel(
            x1'(t) = -(a01 + a21) * x1(t) + a12 * x2(t) + u(t),
            x2'(t) = a21 * x1(t) - a12 * x2(t),
            y(t) = x2(t)
        )
    ),
    # Current result: -delta - beta, delta*beta, -sigma, -c, -b
    # Simpler version: delta + beta, delta*beta, sigma, c, b
    Dict(
        :name => "Goodwin",
        :ode => @ODEmodel(
            x1'(t) = -b * x1(t) + 1 / (c + x4(t)),
            x2'(t) = alpha * x1(t) - beta * x2(t),
            x3'(t) = gama * x2(t) - delta * x3(t),
            x4'(t) = sigma * x4(t) * (gama * x2(t) - delta * x3(t)) / x3(t),
            y(t) = x1(t)
        )
    ),
    # Current result: -M^2, -nu - 432//491*mu, -g, -mu, -b0
    # Simpler version: M^2, mu, g, mu, b0
    Dict(
        :name => "SIRS forced",
        :ode => @ODEmodel(
            s'(t) = mu - mu * s(t) - b0 * (1 + b1 * x1(t)) * i(t) * s(t) + g * r(t),
            i'(t) = b0 * (1 + b1 * x1(t)) * i(t) * s(t) - (nu + mu) * i(t),
            r'(t) = nu * i(t) - (mu + g) * r(t),
            x1'(t) = -M * x2(t),
            x2'(t) = M * x1(t),
            y1(t) = i(t),
            y2(t) = r(t)
        )
    ),
    # Current result: (lm*q - 1)//(lm*d*k^2*u), -1//(c*h*lm*d*k*beta^2), -k*q*beta, -b - 2*h - 2*a, -a, -u, -d, -h
    # Simpler version: (lm*q - 1)//(lm*k^2), c*lm*d*k*beta^2, k*q*beta, b, a, u, d, h (Could be simpler ?)
    Dict(
        :name => "HIV",
        :ode => @ODEmodel(
            x'(t) = lm - d * x(t) - beta * x(t) * v(t),
            y'(t) = beta * x(t) * v(t) - a * y(t),
            v'(t) = k * y(t) - u * v(t),
            w'(t) = c * x(t) * y(t) * w(t) - c * q * y(t) * w(t) - b * w(t),
            z'(t) = c * q * y(t) * w(t) - h * z(t),
            y1(t) = w(t),
            y2(t) = z(t)
        )
    ),
    # The most simplified version of the result I could get by hand is
    # Ninv, b, s, g + a,  (1 - e) * g * (a - s), a * e * g
    # but I have no idea how I did this...
    Dict(
        :name => "SLIQR",
        :ode => @ODEmodel(
            S'(t) = -b * In(t) * S(t) * Ninv - u(t) * S(t) * Ninv,
            L'(t) = b * In(t) * S(t) * Ninv - a * L(t),
            In'(t) = a * L(t) - g * In(t) + s * Q(t),
            Q'(t) = (1 - e) * g * In(t) - s * Q(t),
            y(t) = In(t) * Ninv
        )
    ),
    Dict(
        :name => "St",
        :ode => @ODEmodel(
            S'(t) = r * S(t) - (e + a * W(t)) * S(t) - d * W(t) * S(t) + g * R(t),
            R'(t) = rR * R(t) + (e + a * W(t)) * S(t) - dr * W(t) * R(t) - g * R(t),
            W'(t) = Dd * (T - W(t)),
            y1(t) = S(t) + R(t),
            y2(t) = T
        )
    )
]

# if empty, runs all
to_run = ["St"]

for m in models
    if length(to_run) > 0 && !(m[:name] in to_run)
        continue
    end
    @info "Processing $(m[:name])"
    ode = m[:ode]
    ioeqs = find_ioequations(ode)
    identifiable_functions_raw = extract_identifiable_functions_raw(collect(values(ioeqs)), ode.parameters)

    gens = Array{AbstractAlgebra.Generic.Frac{Nemo.fmpq_mpoly}, 1}(identifiable_functions_raw)
    
    ideal = ParamPunPam.generators_to_saturated_ideal(gens)
    @show ideal
    
    gb = ParamPunPam.paramgb(ideal)
    @show parent(first(gb))
    @show gb
    @info "The coefficients are:"
    for p in gb
        @info collect(coefficients(p))
    end 
end
