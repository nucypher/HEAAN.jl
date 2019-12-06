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


function sub(cipher1::Ciphertext, cipher2::Ciphertext)
    @assert compatible(cipher1, cipher2)
    Ciphertext(
        cipher1.params,
        cipher1.ax - cipher2.ax,
        cipher1.bx - cipher2.bx,
        cipher1.log_cap,
        cipher1.log_precision,
        cipher1.slots)
end


function mul(key::MultiplicationKey, cipher1::Ciphertext, cipher2::Ciphertext)
    # TODO: technically, log_precision may differ?
    @assert compatible(cipher1, cipher2, different_precision=true)
    # TODO: check compatibility with the public key as well

    params = cipher1.params
    plan = rns_plan(params)

    np = cld(2 + cipher1.log_cap + cipher2.log_cap + params.log_polynomial_length + 2, 59)

    ra1 = ntt_rns(to_rns(plan, cipher1.ax, np))
    rb1 = ntt_rns(to_rns(plan, cipher1.bx, np))
    ra2 = ntt_rns(to_rns(plan, cipher2.ax, np))
    rb2 = ntt_rns(to_rns(plan, cipher2.bx, np))

    tp = Polynomial{BinModuloInt{BigInt, cipher1.log_cap}}
    axax = from_rns(tp, ntt_rns(ra1 * ra2, inverse=true))
    bxbx = from_rns(tp, ntt_rns(rb1 * rb2, inverse=true))

    ra1 = ra1 + rb1
    ra2 = ra2 + rb2

    axbx = from_rns(tp, ntt_rns(ra1 * ra2, inverse=true))

    key = key.key

    np = cld(
        cipher1.log_cap + params.log_lo_modulus + params.log_hi_modulus +
            params.log_polynomial_length + 2, 59)

    raa = ntt_rns(to_rns(plan, axax, np))

    #=
    In the paper, we want to find `(P^(-1) * x * y) mod q`,
    where `P = 2^log_hi_modulus`, `Q = 2^log_lo_modulus`, `x` is a polynomial modulo `q <= Q`,
    and `y` is a polynomial modulo `P * Q`.

    So, we perform the multiplication using the full range `q * Q * P`
    (times the polynomial length to prevent overflow during NTT in RNS),
    then drop the high `Q` bits (during the conversion from RNS to big integer),
    then shift by `P`.
    =#

    tp_big = Polynomial{BinModuloInt{BigInt, cipher1.log_cap + params.log_hi_modulus}}
    ax = from_rns(tp_big, ntt_rns(raa * key.rax, inverse=true)) >> params.log_hi_modulus
    bx = from_rns(tp_big, ntt_rns(raa * key.rbx, inverse=true)) >> params.log_hi_modulus

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


function square(mk::MultiplicationKey, cipher::Ciphertext)

    params = cipher.params
    plan = rns_plan(params)

    np = cld(2 + cipher.log_cap * 2 + params.log_polynomial_length + 2, 59)

    ra = ntt_rns(to_rns(plan, cipher.ax, np))
    rb = ntt_rns(to_rns(plan, cipher.bx, np))

    tp = Polynomial{BinModuloInt{BigInt, cipher.log_cap}}
    axax = from_rns(tp, ntt_rns(ra * ra, inverse=true))
    bxbx = from_rns(tp, ntt_rns(rb * rb, inverse=true))
    axbx = from_rns(tp, ntt_rns(ra * rb, inverse=true))
    axbx = axbx + axbx

    key = mk.key

    np = cld(
        cipher.log_cap + params.log_lo_modulus + params.log_hi_modulus +
            params.log_polynomial_length + 2, 59)

    raa = ntt_rns(to_rns(plan, axax, np))

    tp_big = Polynomial{BinModuloInt{BigInt, cipher.log_cap + params.log_hi_modulus}}
    ax = from_rns(tp_big, ntt_rns(raa * key.rax, inverse=true)) >> params.log_hi_modulus
    bx = from_rns(tp_big, ntt_rns(raa * key.rbx, inverse=true)) >> params.log_hi_modulus

    ax = ax + axbx
    bx = bx + bxbx

    Ciphertext(
        params,
        ax,
        bx,
        cipher.log_cap,
        cipher.log_precision * 2,
        cipher.slots)
end


function negate(cipher::Ciphertext)
    Ciphertext(
        cipher.params,
        -cipher.ax,
        -cipher.bx,
        cipher.log_cap,
        cipher.log_precision,
        cipher.slots)
end


function mod_down_by(cipher::Ciphertext, dlog_cap::Int)
    Ciphertext(
        cipher.params,
        mod_down_by(cipher.ax, dlog_cap),
        mod_down_by(cipher.bx, dlog_cap),
        cipher.log_cap - dlog_cap,
        cipher.log_precision,
        cipher.slots)
end


function mod_down_to(cipher::Ciphertext, log_cap::Int)
    Ciphertext(
        cipher.params,
        mod_down_to(cipher.ax, log_cap),
        mod_down_to(cipher.bx, log_cap),
        log_cap,
        cipher.log_precision,
        cipher.slots)
end


function rescale_by(cipher::Ciphertext, dlog_cap::Int)
    Ciphertext(
        cipher.params,
        cipher.ax >> dlog_cap,
        cipher.bx >> dlog_cap,
        cipher.log_cap - dlog_cap,
        cipher.log_precision - dlog_cap,
        cipher.slots)
end


function add_const(cipher::Ciphertext, cnst::Float64)
    tp = BinModuloInt{BigInt, cipher.log_cap}
    cnst_big = float_to_integer(tp, cnst, cipher.log_precision)
    Ciphertext(
        cipher.params,
        cipher.ax,
        cipher.bx + cnst_big,
        cipher.log_cap,
        cipher.log_precision,
        cipher.slots)
end


function add_const(cipher::Ciphertext, cnst::Complex{Float64})
    negate(imul(add_const(imul(add_const(cipher, real(cnst))), -imag(cnst))))
end


function mul_by_const(cipher::Ciphertext, cnst::Float64, log_precision::Int)
    tp = BinModuloInt{BigInt, cipher.log_cap}
    cnst_big = float_to_integer(tp, cnst, log_precision)
    Ciphertext(
        cipher.params,
        cipher.ax * cnst_big,
        cipher.bx * cnst_big,
        cipher.log_cap,
        cipher.log_precision + log_precision,
        cipher.slots)
end


function mul_by_const(cipher::Ciphertext, cnst::Complex{Float64}, log_precision::Int)
    add(
        mul_by_const(cipher, real(cnst), log_precision),
        imul(mul_by_const(cipher, imag(cnst), log_precision)))
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


function mul_by_const_vec(cipher::Ciphertext, v::Array{Complex{Float64}, 1}, log_precision::Int)
    p = encode(cipher.params, v, log_precision, cipher.log_cap)
    mul_by_plaintext(cipher, p)
end


function mul_by_plaintext(cipher::Ciphertext, p::Plaintext)
    params = cipher.params
    bnd = maximum(num_bits.(p.polynomial.coeffs))
    np = cld(cipher.log_cap + bnd + params.log_polynomial_length + 2, 59)

    plan = rns_plan(params)
    rpoly = ntt_rns(to_rns(plan, p.polynomial, np))

    # TODO: use mul_by_rns() here
    Ciphertext(
        cipher.params,
        mult(cipher.ax, rpoly, np),
        mult(cipher.bx, rpoly, np),
        cipher.log_cap,
        cipher.log_precision + p.log_precision,
        cipher.slots)
end


function Base.circshift(rk::LeftRotationKey, cipher::Ciphertext, shift::Integer)
    params = cipher.params
    plan = rns_plan(params)

    # TODO: handle positive shifts too (mod N?)
    shift = -shift

    @assert rk.shift == shift

    axrot = left_rotate(cipher.ax, shift)
    bxrot = left_rotate(cipher.bx, shift)

    np = cld(cipher.log_cap + params.log_lo_modulus +
        params.log_hi_modulus + params.log_polynomial_length + 2, 59)

    rarot = ntt_rns(to_rns(plan, axrot, np))

    tp_big = Polynomial{BinModuloInt{BigInt, cipher.log_cap + params.log_hi_modulus}}
    ax = from_rns(tp_big, ntt_rns(rarot * rk.key.rax, inverse=true)) >> params.log_hi_modulus
    bx = from_rns(tp_big, ntt_rns(rarot * rk.key.rbx, inverse=true)) >> params.log_hi_modulus
    bx = bx + bxrot

    Ciphertext(
        params,
        ax,
        bx,
        cipher.log_cap,
        cipher.log_precision,
        cipher.slots)
end


function Base.conj(ck::ConjugationKey, cipher::Ciphertext)
    params = cipher.params
    plan = rns_plan(params)

    axconj = conjugate(cipher.ax)
    bxconj = conjugate(cipher.bx)

    np = cld(cipher.log_cap + params.log_lo_modulus +
        params.log_hi_modulus + params.log_polynomial_length + 2, 59)
    raconj = ntt_rns(to_rns(plan, axconj, np))

    tp_big = Polynomial{BinModuloInt{BigInt, cipher.log_cap + params.log_hi_modulus}}
    ax = from_rns(tp_big, ntt_rns(raconj * ck.key.rax, inverse=true)) >> params.log_hi_modulus
    bx = from_rns(tp_big, ntt_rns(raconj * ck.key.rbx, inverse=true)) >> params.log_hi_modulus
    bx = bx + bxconj

    Ciphertext(
        params,
        ax,
        bx,
        cipher.log_cap,
        cipher.log_precision,
        cipher.slots)
end
