function add(cipher1::Ciphertext, cipher2::Ciphertext)
    @assert compatible(cipher1, cipher2)
    Ciphertext(
        cipher1.params,
        cipher1.ax + cipher2.ax,
        cipher1.bx + cipher2.bx,
        cipher1.log_cap,
        cipher1.log_precision,
        cipher1.slots)
end


function mul(pk::PublicKeySet, cipher1::Ciphertext, cipher2::Ciphertext)
    # TODO: technically, log_precision may differ?
    @assert compatible(cipher1, cipher2)
    # TODO: check compatibility with the public key as well

    params = cipher1.params
    plan = rns_plan(params)

    np = cld(2 + cipher1.log_cap + cipher2.log_cap + params.log_polynomial_length + 2, 59)

    ra1 = ntt_rns(to_rns(plan, cipher1.ax, np))
    rb1 = ntt_rns(to_rns(plan, cipher1.bx, np))
    ra2 = ntt_rns(to_rns(plan, cipher2.ax, np))
    rb2 = ntt_rns(to_rns(plan, cipher2.bx, np))

    tp = Polynomial{BinModuloInt{BigInt, cipher1.log_cap}}
    axax = from_rns(tp, ntt_rns(ra1 * ra2, inverse=true), np)
    bxbx = from_rns(tp, ntt_rns(rb1 * rb2, inverse=true), np)

    ra1 = ra1 + rb1
    ra2 = ra2 + rb2

    axbx = from_rns(tp, ntt_rns(ra1 * ra2, inverse=true), np)

    key = pk.mul_key.key

    np = cld(
        cipher1.log_cap + params.log_lo_modulus + params.log_hi_modulus +
            params.log_polynomial_length + 2, 59)

    raa = ntt_rns(to_rns(plan, axax, np))

    tp_big = Polynomial{BinModuloInt{BigInt, cipher1.log_cap + params.log_lo_modulus}}
    ax = from_rns(tp_big, ntt_rns(raa * key.rax, inverse=true), np) >> params.log_lo_modulus
    bx = from_rns(tp_big, ntt_rns(raa * key.rbx, inverse=true), np) >> params.log_lo_modulus

    ax = ax + axbx - bxbx - axax
    bx = bx + bxbx

    Ciphertext(
        params,
        ax,
        bx,
        cipher1.log_cap,
        cipher1.log_precision + cipher2.log_precision,
        cipher1.slots)
end


# TODO: do we even need this as a separate function?
function mul_by_monomial(x::Polynomial, pwr::Integer)
    shift_polynomial(x, pwr)
end


function imul(cipher::Ciphertext)
    params = cipher.params
    Ciphertext(
        params,
        mul_by_monomial(cipher.ax, 2^(params.log_polynomial_length - 1)),
        mul_by_monomial(cipher.bx, 2^(params.log_polynomial_length - 1)),
        cipher.log_cap,
        cipher.log_precision,
        cipher.slots)
end
