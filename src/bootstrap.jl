struct BootContext

    rpvec :: Array{RNSPolynomialTransformed, 1}
    rpvec_inv :: Array{RNSPolynomialTransformed, 1}
    rp1 :: RNSPolynomialTransformed
    rp2 :: RNSPolynomialTransformed

    log_precision :: Int
    log_slots :: Int

    function BootContext(params::Params, log_slots::Int, log_precision::Int)

        e_plan = embedding_plan(params)
        r_plan = rns_plan(params)

        slots = 1 << log_slots
        dslots = slots << 1
        logk = log_slots >> 1

        k = 1 << logk
        gap = 2^(params.log_polynomial_length - 1) >> log_slots

        log_plen = params.log_polynomial_length
        plen = 2^log_plen

        rpvec = Array{RNSPolynomialTransformed}(undef, slots)
        rpvec_inv = Array{RNSPolynomialTransformed}(undef, slots)

        tp = BinModuloInt{BigInt, log_precision + 1}

        c = 0.25 / pi

        if log_slots < log_plen - 1
            dgap = gap >> 1
            for pos in 0:slots-1
                ki = pos ÷ k * k
                i = collect(0:slots-1)
                deg = mod.(.-(e_plan.rotation_group[(i .+ pos) .% slots .+ 1]) .* i .* gap, 2 * plen)
                pvals = vcat(e_plan.root_powers[deg .+ 1], e_plan.root_powers[deg .+ 1] .* im)
                pvals = circshift(pvals, ki)
                pvals = unembed(e_plan, pvals)
                pvals_real = vcat(real.(pvals), imag.(pvals))

                pvec = zeros(tp, plen)
                pvec[1:dgap:end] = float_to_integer.(tp, pvals_real, log_precision)

                rpvec[pos + 1] = to_rns_transformed(
                    r_plan, Polynomial(pvec, negacyclic_modulus),
                    params.log_lo_modulus + 2 * log_plen)
            end

            pvals = vcat([0.0 for i in 1:slots], [-c*im for i in 1:slots])
            pvals = unembed(e_plan, pvals)
            pvals_real = vcat(real.(pvals), imag.(pvals))

            pvec = zeros(tp, plen)
            pvec[1:dgap:end] = float_to_integer.(tp, pvals_real, log_precision)

            rp1 = to_rns_transformed(
                r_plan, Polynomial(pvec, negacyclic_modulus),
                params.log_lo_modulus + 2 * log_plen)


            pvals = vcat([complex(c) for i in 1:slots], [0.0 for i in 1:slots])
            pvals = unembed(e_plan, pvals)
            pvals_real = vcat(real.(pvals), imag.(pvals))

            pvec = zeros(tp, plen)
            pvec[1:dgap:end] = float_to_integer.(tp, pvals_real, log_precision)

            rp2 = to_rns_transformed(
                r_plan, Polynomial(pvec, negacyclic_modulus),
                params.log_lo_modulus + 2 * log_plen)
            for i in 0:plen-1
                pvec[i+1] = zero(tp)
            end
        else
            # TODO: (issue #1) need to test this branch
            for pos in 0:slots-1
                ki = pos ÷ k * k

                i = collect(0:slots-1)
                deg = mod.(.-(e_plan.rotation_group[(i .+ pos) .% slots .+ 1]) .* i .* gap, plen * 2)
                pvals = e_plan.root_powers[deg.+1]

                pvals = circshift(pvals[1:slots], ki)
                pvals = unembed(e_plan, pvals[1:slots])
                pvals_real = vcat(real.(pvals), imag.(pvals))

                pvec = zeros(tp, plen)
                pvec[1:gap:end] = float_to_integer.(tp, pvals_real, log_precision)

                rpvec[pos+1] = to_rns_transformed(
                    r_plan, Polynomial(pvec, negacyclic_modulus),
                    params.log_lo_modulus + 2 * log_plen)
            end

            # These will be unused
            rp1 = RNSPolynomialTransformed(r_plan, zeros(UInt64, plen, 1), 0, 0, false)
            rp2 = RNSPolynomialTransformed(r_plan, zeros(UInt64, plen, 1), 0, 0, false)
        end

        for pos in 0:slots-1
            ki = pos ÷ k * k

            i = collect(0:slots-1)
            deg = (e_plan.rotation_group[i.+1] .* ((i .+ pos) .% slots) .* gap) .% (plen * 2)
            pvals = e_plan.root_powers[deg.+1]
            pvals = circshift(pvals[1:slots], ki)
            pvals = unembed(e_plan, pvals[1:slots])

            pvals_real = vcat(real.(pvals), imag.(pvals))

            pvec = zeros(tp, plen)
            pvec[1:gap:end] = float_to_integer.(tp, pvals_real, log_precision)

            rpvec_inv[pos+1] = to_rns_transformed(
                r_plan, Polynomial(pvec, negacyclic_modulus),
                params.log_lo_modulus + 2 * log_plen)
        end

        new(rpvec, rpvec_inv, rp1, rp2, log_precision, log_slots)
    end
end


"""
    BootstrapKey(
        rng::AbstractRNG, secret_key::SecretKey,
        enc_key::EncryptionKey, mul_key::MultiplicationKey, conj_key::ConjugationKey,
        log_slots::Int, log_precision::Int)

A public key used for bootstrapping.
Can only be applied to encrypted vector of size `2^log_slots`.
Needs the [`SecretKey`](@ref), an [`EncryptionKey`](@ref),
a [`MultiplicationKey`](@ref) and a [`ConjugationKey`](@ref).
"""
struct BootstrapKey

    enc_key :: EncryptionKey
    mul_key :: MultiplicationKey
    conj_key :: ConjugationKey
    rot_keys :: Dict{Int, LeftRotationKey}
    log_slots :: Int
    bc :: BootContext

    function BootstrapKey(
            rng::AbstractRNG, secret_key::SecretKey,
            enc_key::EncryptionKey, mul_key::MultiplicationKey, conj_key::ConjugationKey,
            log_slots::Int, log_precision::Int)

        params = secret_key.params
        log_plen = params.log_polynomial_length

        loglh = log_slots ÷ 2
        k = 1 << loglh
        m = 1 << (log_slots - loglh)

        shifts = Set(vcat(
            1 .<< collect(0:log_plen-2),
            collect(1:k-1),
            collect(1:m-1) .* k
            ))

        rot_keys = Dict(shift => LeftRotationKey(rng, secret_key, shift) for shift in shifts)

        bc = BootContext(params, log_slots, log_precision)

        new(enc_key, mul_key, conj_key, rot_keys, log_slots, bc)
    end
end


left_rotate(bk::BootstrapKey, c::Ciphertext, shift::Int) = circshift(bk.rot_keys[shift], c, -shift)


function _coeff_to_slot(bk::BootstrapKey, bc::BootContext, cipher::Ciphertext, inverse::Bool)

    # `inverse = false` - coeffToSlot in the paper
    # `inverse = true`  - slotToCoeff in the paper

    @assert 2^bc.log_slots == cipher.slots

    rpvec = inverse ? bc.rpvec_inv : bc.rpvec

    slots = cipher.slots
    log_slots = bc.log_slots
    logk = log_slots ÷ 2
    k = 1 << logk

    rotvec = Array{Ciphertext}(undef, k)
    rotvec[1] = cipher
    for j in 2:k
        rotvec[j] = left_rotate(bk, cipher, j-1)
    end

    tmpvec = Array{Ciphertext}(undef, k)
    for j in 1:k
        tmpvec[j] = mul_by_rns(rotvec[j], rpvec[j], bc.log_precision)
    end
    for j in 2:k
        tmpvec[1] = add(tmpvec[1], tmpvec[j])
    end

    cipher = tmpvec[1]

    for ki in k:k:slots-1
        for j in 1:k
            tmpvec[j] = mul_by_rns(rotvec[j], rpvec[j+ki], bc.log_precision)
        end
        for j in 2:k
            tmpvec[1] = add(tmpvec[1], tmpvec[j])
        end
        tmpvec[1] = left_rotate(bk, tmpvec[1], ki)
        cipher = add(cipher, tmpvec[1])
    end

    rescale_by(cipher, bc.log_precision)
end


coeff_to_slot(bk, bc, cipher) = _coeff_to_slot(bk, bc, cipher, false)
slot_to_coeff(bk, bc, cipher) = _coeff_to_slot(bk, bc, cipher, true)


function exp2pi(mk::MultiplicationKey, cipher::Ciphertext, log_precision::Int)
    Pi = Float64(pi)

    # Cipher: log_precision, log_cap = (p, q); log_precision = l

    cipher2 = square(mk, cipher) # (2p, q)
    cipher2 = rescale_by(cipher2, log_precision) # (2p - l, q - l)

    cipher4 = square(mk, cipher2) # (4p - 2l, q - l)
    cipher4 = rescale_by(cipher4, log_precision) # (4p - 3l, q - 2l)
    c = 1/(2 * Pi)
    cipher01 = add_const(cipher, c) # (p, q)

    c = 2*Pi
    cipher01 = mul_by_const(cipher01, c, log_precision) # (p + l, q)
    cipher01 = rescale_by(cipher01, log_precision) # (p, q - l)

    c = 3/(2*Pi)
    cipher23 = add_const(cipher, c) # (p, q)

    c = 4*Pi*Pi*Pi/3
    cipher23 = mul_by_const(cipher23, c, log_precision) # (p + l, q)
    cipher23 = rescale_by(cipher23, log_precision) # (p, q - l)

    cipher23 = mul(mk, cipher23, cipher2) # (p, q - l) * (2p - l, q - l) = (3p - l, q - l)
    cipher23 = rescale_by(cipher23, log_precision) # (3p - 2l, q - 2l)

    # TODO: (issue #1) how justified is mod_down_to() here? In the original these two
    # ciphertexts are just added without regard to different log_caps, and the log_cap of
    # cipher23 is used for the result (which is smaller than that of cipher01).
    # So we have an inconsistent Ciphertext in the output, but somehow it works out...
    # (3p - 2l, q - 2l) + (p, q)
    cipher23 = add(cipher23, mod_down_to(cipher01, cipher23.log_cap))

    c = 5/(2*Pi)
    cipher45 = add_const(cipher, c) # (p, q)

    c = 4*Pi*Pi*Pi*Pi*Pi/15
    cipher45 = mul_by_const(cipher45, c, log_precision)
    cipher45 = rescale_by(cipher45, log_precision) # cipher45.log_cap : log_cap - log_precision

    c = 7/(2*Pi)
    cipher = add_const(cipher, c) # cipher.log_cap : log_cap

    c = 8*Pi*Pi*Pi*Pi*Pi*Pi*Pi/315
    cipher = mul_by_const(cipher, c, log_precision)
    cipher = rescale_by(cipher, log_precision) # cipher.log_cap : log_cap - log_precision

    cipher = mul(mk, cipher, cipher2)
    cipher = rescale_by(cipher, log_precision) # cipher.log_cap : log_cap - 2log_precision

    cipher45 = mod_down_by(cipher45, log_precision) # cipher45.log_cap : log_cap - 2log_precision
    cipher = add(cipher, cipher45) # cipher.log_cap : log_cap - 2log_precision

    cipher = mul(mk, cipher, cipher4)
    cipher = rescale_by(cipher, log_precision) # cipher.log_cap : log_cap - 3log_precision

    cipher23 = mod_down_by(cipher23, log_precision)
    cipher = add(cipher, cipher23) # cipher.log_cap : log_cap - 3log_precision

    cipher
end


function eval_exp(bk::BootstrapKey, bc::BootContext, cipher::Ciphertext, log_t::Int, log_i::Int=4)

    @assert 2^bc.log_slots == cipher.slots
    slots = cipher.slots
    log_slots = bc.log_slots

    mk = bk.mul_key
    ck = bk.conj_key

    log_plen = cipher.params.log_polynomial_length

    if log_slots < log_plen - 1
        tmp = conj(ck, cipher)
        cipher = sub(cipher, tmp)
        cipher = div_by_po2(cipher, log_t + 1) # bitDown: log_t + 1
        cipher = exp2pi(mk, cipher, bc.log_precision) # bitDown: log_t + 1 + 3(log_cap + log_i)
        for i in 0:log_i+log_t-1
            cipher = square(mk, cipher)
            cipher = rescale_by(cipher, bc.log_precision)
        end
        tmp = conj(ck, cipher)
        cipher = sub(cipher, tmp)
        tmp = mul_by_rns(cipher, bc.rp1, bc.log_precision)
        tmprot = circshift(bk.rot_keys[slots], tmp, -slots)
        tmp = add(tmp, tmprot)
        cipher = mul_by_rns(cipher, bc.rp2, bc.log_precision)
        tmprot = circshift(bk.rot_keys[slots], cipher, -slots)
        cipher = add(cipher, tmprot)
        cipher = add(cipher, tmp)
    else # TODO: (issue #1) check this branch
        tmp = conj(ck, cipher)
        c2 = sub(cipher, tmp)
        cipher = add(cipher, tmp)
        cipher = imul(cipher)
        cipher = div_by_po2(cipher, log_t + 1) # cipher bitDown: log_t + 1
        c2 = rescale_by(c2, log_t + 1) # c2 bitDown: log_t + 1
        cipher = exp2pi(mk, cipher, bc.log_precision) # cipher bitDown: log_t + 1 + 3(log_cap + log_i)
        c2 = exp2pi(mk, c2, bc.log_precision) # c2 bitDown: log_t + 1 + 3(log_cap + log_i)
        for i in 0:log_i+log_t-1
            c2 = square(mk, c2)
            cipher = square(mk, cipher)
            c2 = rescale_by(c2, bc.log_precision)
            cipher = rescale_by(cipher, bc.log_precision)
        end
        tmp = conj(ck, c2)
        c2 = sub(c2, tmp)
        tmp = conj(ck, cipher)
        cipher = sub(cipher, tmp)
        cipher = imul(cipher)
        cipher = sub(c2, cipher)
        c = 0.25/pi
        cipher = mul_by_const(cipher, c, bc.log_precision)
    end
    rescale_by(cipher, bc.log_precision + log_i)
end


#=
The theory behind this is explained in
"Improved Bootstrapping for Approximate Homomorphic Encryption", Cheon et al (2018),
Section 2.2.
=#

"""
    bootstrap(bk::BootstrapKey, cipher::Ciphertext, log_t::Int=4)

Bootstrap a ciphertext, increasing the difference between its `log_cap` and `log_precision`.
Needs a [`BootstrapKey`](@ref) object.
"""
function bootstrap(bk::BootstrapKey, cipher::Ciphertext, log_t::Int=4)

    @assert 2^bk.log_slots == cipher.slots
    params = cipher.params

    log_plen = params.log_polynomial_length
    orig_log_precision = cipher.log_precision

    # TODO: (issue #11) check that bk.bc.log_precision >= cipher.log_precision?
    cipher = Ciphertext(
        params,
        mod_up_to(cipher.ax, params.log_lo_modulus),
        mod_up_to(cipher.bx, params.log_lo_modulus),
        params.log_lo_modulus,
        bk.bc.log_precision,
        cipher.slots)

    for i in bk.log_slots:log_plen-2
        rot = circshift(bk.rot_keys[1 << i], cipher, -(1 << i))
        cipher = add(cipher, rot)
    end

    cipher = div_by_po2(cipher, log_plen - 1)
    cipher = coeff_to_slot(bk, bk.bc, cipher)
    cipher = eval_exp(bk, bk.bc, cipher, log_t)
    cipher = slot_to_coeff(bk, bk.bc, cipher)

    Ciphertext(
        params,
        cipher.ax,
        cipher.bx,
        cipher.log_cap,
        orig_log_precision,
        cipher.slots)
end
