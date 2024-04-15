
mutable struct GroebnerState{Blackbox, FF, PolyFF, PolyFracQQ, OrderingGb}
    # original polynomials over Q(a)
    blackbox::Blackbox
    shape::Vector{Vector{PolyFF}}
    # total degrees of the parameters of the Groebner basis
    param_degrees::Vector{Vector{Tuple{Int, Int}}}
    # exponents of the parameters of the Groebner basis
    param_exponents::Vector{Vector{Tuple{PolyFF, PolyFF}}}
    param_coeffs_crt::Vector{Vector{Tuple{Vector{BigInt}, Vector{BigInt}}}}
    field_to_param_exponents::Dict{FF, Vector{Vector{Tuple{PolyFF, PolyFF}}}}
    # fully reconstructed basis
    param_basis::Vector{PolyFracQQ}
    # GB computation helpers
    gb_context::Any
    gb_ordering::OrderingGb

    function GroebnerState(
        blackbox::Blackbox,
        ord::Ord
    ) where {Blackbox <: AbstractBlackboxIdeal, Ord}
        Rx = parent(blackbox)
        Ra = parent_params(blackbox)
        params = gens(Ra)
        polyvars = gens(Rx)
        K = base_ring(Ra)
        @debug "Given $(length(blackbox)) functions in $K($(join(repr.(params),", ")))[$(join(repr.(polyvars),", "))]"
        PolyFF = Nemo.fpMPolyRingElem
        PolyFracQQ = Any
        FF = Any
        new{Blackbox, FF, PolyFF, PolyFracQQ, typeof(ord)}(
            blackbox,
            Vector{Vector{Nemo.fpMPolyRingElem}}(),
            Vector{Vector{Tuple{Int, Int}}}(),
            Vector{Vector{Tuple{PolyFF, PolyFF}}}(),
            Vector{Vector{Tuple{Vector{BigInt}, Vector{BigInt}}}}(),
            Dict{FF, Vector{Vector{Tuple{PolyFF, PolyFF}}}}(),
            Vector{PolyFracQQ}(),
            nothing,
            ord
        )
    end
end

function randluckyspecpoint(state::GroebnerState, K)
    # TODO: this is not correct!
    Ra = parent_params(state.blackbox)
    [rand(K) for _ in 1:nvars(Ra)]
end

function reconstruct_crt!(state, modular)
    field_to_param_exponents = state.field_to_param_exponents
    char = UInt64(characteristic(modular.finite_field))
    if length(field_to_param_exponents) == 1
        shape = state.shape
        param_coeffs_crt =
            Vector{Vector{Tuple{Vector{BigInt},Vector{BigInt}}}}(undef, length(shape))
        for i in 1:length(param_coeffs_crt)
            param_coeffs_crt[i] =
                Vector{Tuple{Vector{BigInt},Vector{BigInt}}}(undef, length(shape[i]))
            for j in 1:length(param_coeffs_crt[i])
                P, Q = state.param_exponents[i][j]
                Pcoeffs = map(c -> BigInt(data(c)), collect(coefficients(P)))
                Qcoeffs = map(c -> BigInt(data(c)), collect(coefficients(Q)))
                param_coeffs_crt[i][j] = (Pcoeffs, Qcoeffs)
            end
        end
        @debug "CRT-Reconstructed coefficients" param_coeffs_crt
        state.param_coeffs_crt = param_coeffs_crt
        modular.modulo *= char
        return nothing
    end
    buf, n1, n2, M, bigch, invm1, invm2 =
        BigInt(), BigInt(), BigInt(), BigInt(), BigInt(), BigInt(), BigInt()
    # Base.GMP.MPZ.set_ui!(bigch, char)
    # Base.GMP.MPZ.mul_ui!(M, modular.modulo, char)
    # Base.GMP.MPZ.gcdext!(buf, invm1, invm2, modular.modulo, bigch)
    Groebner.crt_precompute!(M, n1, n2, invm1, modular.modulo, invm2, char)
    param_exponents = state.param_exponents
    param_coeffs_crt = state.param_coeffs_crt
    for i in 1:length(param_exponents)
        for j in 1:length(param_exponents[i])
            @assert length(param_exponents[i]) == length(param_coeffs_crt[i])
            for k in 1:length(param_exponents[i][j][1])
                ca = param_coeffs_crt[i][j][1][k]
                cf = UInt64(data(coeff(param_exponents[i][j][1], k)))
                Groebner.crt!(M, buf, n1, n2, ca, invm1, cf, invm2)
                Base.GMP.MPZ.set!(param_coeffs_crt[i][j][1][k], buf)
            end
            for k in 1:length(param_exponents[i][j][2])
                ca = param_coeffs_crt[i][j][2][k]
                cf = UInt(data(coeff(param_exponents[i][j][2], k)))
                Groebner.crt!(M, buf, n1, n2, ca, invm1, cf, invm2)
                Base.GMP.MPZ.set!(param_coeffs_crt[i][j][2][k], buf)
            end
        end
    end
    @debug "CRT-Reconstructed coefficients" param_coeffs_crt
    modular.modulo *= char
    nothing
end

"""
    rational_reconstruct_polynomial(ring, poly_ff)

Reconstructs each coefficient of the given polynomial into rationals.
Resulting polynomial lives in the given `ring`.

Returns (success, poly_qq), where `success` is `true` if the reconstruction was
successful and `false`, otherwise.
"""
function rational_reconstruct_polynomial(ring, poly_ff)
    finite_field = base_ring(parent(poly_ff))
    modulo = BigInt(characteristic(finite_field))
    bnd = Groebner.ratrec_reconstruction_bound(modulo)
    buf, buf1 = BigInt(), BigInt()
    buf2, buf3 = BigInt(), BigInt()
    u1, u2 = BigInt(), BigInt()
    u3, v1 = BigInt(), BigInt()
    v2, v3 = BigInt(), BigInt()
    cfs_ff = collect(coefficients(poly_ff))
    cfs_qq = map(_ -> Rational{BigInt}(0), cfs_ff)
    cfs_rec = map(_ -> Nemo.QQ(0), cfs_ff)
    evs = collect(exponent_vectors(poly_ff))
    success = true
    for i in 1:length(cfs_ff)
        cz = BigInt(data(cfs_ff[i]))
        cq = cfs_qq[i]
        num, den = numerator(cq), denominator(cq)
        success_ = Groebner.ratrec!(
            num,
            den,
            bnd,
            buf,
            buf1,
            buf2,
            buf3,
            u1,
            u2,
            u3,
            v1,
            v2,
            v3,
            cz,
            modulo
        )
        cfs_rec[i] = Nemo.QQ(cfs_qq[i])
        success = success && success_
    end
    poly_qq = ring(cfs_rec, evs)
    return success, poly_qq
end

function reconstruct_rational!(state, modular)
    blackbox = state.blackbox
    Rorig = parent(blackbox)
    Rparam = parent_params(blackbox)
    Rorig_frac, _ = polynomial_ring(
        Nemo.fraction_field(Rparam),
        symbols(Rorig),
        internal_ordering=Nemo.internal_ordering(Rorig)
    )
    Rparam_frac = base_ring(Rorig_frac)
    polysreconstructed = Vector{elem_type(Rorig_frac)}(undef, length(state.shape))
    modulo = modular.modulo
    bnd = Groebner.ratrec_reconstruction_bound(modulo)
    buf, buf1 = BigInt(), BigInt()
    buf2, buf3 = BigInt(), BigInt()
    u1, u2 = BigInt(), BigInt()
    u3, v1 = BigInt(), BigInt()
    v2, v3 = BigInt(), BigInt()
    param_coeffs_crt = state.param_coeffs_crt
    param_exponents = state.param_exponents
    @debug "Reconstruction" modulo bnd
    for i in 1:length(param_coeffs_crt)
        coeffsrec = Vector{elem_type(Rparam_frac)}(undef, length(state.shape[i]))
        # skip reconstrction of the first coefficient, it is equal to one in the
        # reduced basis
        for j in 1:length(param_coeffs_crt[i])
            rec_coeffs = Vector{Rational{BigInt}}(undef, length(param_coeffs_crt[i][j][1]))
            for k in 1:length(param_coeffs_crt[i][j][1])
                cz = param_coeffs_crt[i][j][1][k]
                cq = Rational{BigInt}(0)
                num, den = numerator(cq), denominator(cq)
                success = Groebner.ratrec!(
                    num,
                    den,
                    bnd,
                    buf,
                    buf1,
                    buf2,
                    buf3,
                    u1,
                    u2,
                    u3,
                    v1,
                    v2,
                    v3,
                    cz,
                    modulo
                )
                rec_coeffs[k] = cq
                !success && return false
            end
            exponent_vecs = collect(exponent_vectors(param_exponents[i][j][1]))
            param_num = Rparam(rec_coeffs, exponent_vecs)
            rec_coeffs = Vector{Rational{BigInt}}(undef, length(param_coeffs_crt[i][j][2]))
            for k in 1:length(param_coeffs_crt[i][j][2])
                cz = param_coeffs_crt[i][j][2][k]
                cq = Rational{BigInt}(0)
                num, den = numerator(cq), denominator(cq)
                success = Groebner.ratrec!(
                    num,
                    den,
                    bnd,
                    buf,
                    buf1,
                    buf2,
                    buf3,
                    u1,
                    u2,
                    u3,
                    v1,
                    v2,
                    v3,
                    cz,
                    modulo
                )
                rec_coeffs[k] = cq
                !success && return false
            end
            exponent_vecs = collect(exponent_vectors(param_exponents[i][j][2]))
            param_den = Rparam(rec_coeffs, exponent_vecs)
            coeffsrec[j] = param_num // param_den
        end
        @debug "QQ-Reconstructed coefficients" coeffsrec
        polysreconstructed[i] =
            Rorig_frac(coeffsrec, map(e -> exponent_vector(e, 1), state.shape[i]))
    end
    state.param_basis = polysreconstructed
    true
end

function assess_correctness_mod_p(state, modular)
    # NOTE: in this state, ideally, this should always pass, since the modulo
    # we use to check correctness is the same as the one used to compute the
    # basis
    find_next_lucky_prime!(modular)
    finite_field = modular.finite_field
    reduce_mod_p!(state.blackbox, finite_field)
    point = randluckyspecpoint(state, finite_field)
    @debug "Checking correctness at $point in $finite_field"
    generators_zp = specialize_mod_p(state.blackbox, point)
    R_zp = parent(first(generators_zp))
    basis_specialized_coeffs = map(
        f -> map(
            c ->
                evaluate(map_coefficients(cc -> finite_field(cc), numerator(c)), point) // evaluate(
                    map_coefficients(cc -> finite_field(cc), denominator(c)),
                    point
                ),
            collect(coefficients(f))
        ),
        state.param_basis
    )
    param_basis_specialized = map(
        i -> R_zp(
            basis_specialized_coeffs[i],
            collect(exponent_vectors(state.param_basis[i]))
        ),
        1:length(basis_specialized_coeffs)
    )
    @debug "Evaluated basis" param_basis_specialized
    @debug "Evaluated generators" generators_zp
    if !isgroebner(param_basis_specialized, ordering=state.gb_ordering)
        @debug "Not a basis!"
        return false
    end
    inclusion = normalform(param_basis_specialized, generators_zp, ordering=state.gb_ordering)
    @debug "Inclusion in correctness assessment" inclusion
    all(iszero, inclusion)
end
