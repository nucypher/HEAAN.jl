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
    @assert compatible(cipher1, cipher2, different_precision=true)
    # TODO: (issue #11) check compatibility with the public key as well

    params = cipher1.params
    plan = rns_plan(params)

    log_cap = cipher1.log_cap

    ra1 = to_rns_transformed(plan, cipher1.ax, log_cap)
    rb1 = to_rns_transformed(plan, cipher1.bx, log_cap)
    ra2 = to_rns_transformed(plan, cipher2.ax, log_cap)
    rb2 = to_rns_transformed(plan, cipher2.bx, log_cap)

    axax = from_rns_transformed(ra1 * ra2, log_cap)
    bxbx = from_rns_transformed(rb1 * rb2, log_cap)

    ra1 = ra1 + rb1
    ra2 = ra2 + rb2

    axbx = from_rns_transformed(ra1 * ra2, log_cap)

    key = key.key

    raa = to_rns_transformed(plan, axax, params.log_lo_modulus + params.log_hi_modulus)

    #=
    In the paper, we want to find `(P^(-1) * x * y) mod q`,
    where `P = 2^log_hi_modulus`, `Q = 2^log_lo_modulus`, `x` is a polynomial modulo `q <= Q`,
    and `y` is a polynomial modulo `P * Q`.

    So, we perform the multiplication using the full range `q * Q * P`
    (times the polynomial length to prevent overflow during NTT in RNS),
    then drop the high `Q` bits (during the conversion from RNS to big integer),
    then shift by `P`.
    =#

    log_modulus = cipher1.log_cap + params.log_hi_modulus
    ax = from_rns_transformed(raa * key.rax, log_modulus) >> params.log_hi_modulus
    bx = from_rns_transformed(raa * key.rbx, log_modulus) >> params.log_hi_modulus

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
    log_cap = cipher.log_cap

    ra = to_rns_transformed(plan, cipher.ax, log_cap)
    rb = to_rns_transformed(plan, cipher.bx, log_cap)

    axax = from_rns_transformed(ra * ra, log_cap)
    bxbx = from_rns_transformed(rb * rb, log_cap)
    axbx = from_rns_transformed(ra * rb, log_cap)
    axbx = axbx + axbx

    key = mk.key

    raa = to_rns_transformed(plan, axax, params.log_lo_modulus + params.log_hi_modulus)

    log_modulus = cipher.log_cap + params.log_hi_modulus
    ax = from_rns_transformed(raa * key.rax, log_modulus) >> params.log_hi_modulus
    bx = from_rns_transformed(raa * key.rbx, log_modulus) >> params.log_hi_modulus

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


function div_by_po2(cipher::Ciphertext, bits::Int)
    Ciphertext(
        cipher.params,
        cipher.ax >> bits,
        cipher.bx >> bits,
        cipher.log_cap - bits,
        cipher.log_precision,
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
    shift = 1 << (params.log_polynomial_length - 1)
    Ciphertext(
        params,
        shift_polynomial(cipher.ax, shift),
        shift_polynomial(cipher.bx, shift),
        cipher.log_cap,
        cipher.log_precision,
        cipher.slots)
end


function mul_by_const_vec(cipher::Ciphertext, v::Array{Complex{Float64}, 1}, log_precision::Int)
    p = encode(cipher.params, v, log_precision, cipher.log_cap, true)
    mul_by_plaintext(cipher, p)
end


function mul_by_plaintext(cipher::Ciphertext, p::Plaintext)
    @assert cipher.params == p.params
    @assert cipher.slots == p.slots
    # `cipher` and `p` can have different `log_precision` (since it is multiplication),
    # but also different `log_cap` (since for the plaintext we can increase `log_cap` safely).

    params = cipher.params

    plan = rns_plan(params)
    rpoly = to_rns_transformed(plan, p.polynomial, cipher.log_cap)

    mul_by_rns(cipher, rpoly, p.log_precision)
end


function Base.circshift(rk::LeftRotationKey, cipher::Ciphertext, shift::Integer)
    params = cipher.params
    plan = rns_plan(params)

    # TODO: (issue #12) handle positive shifts too (mod N?)
    shift = -shift

    @assert rk.shift == shift

    axrot = left_rotate(cipher.ax, shift)
    bxrot = left_rotate(cipher.bx, shift)

    rarot = to_rns_transformed(plan, axrot, params.log_lo_modulus + params.log_hi_modulus)

    log_modulus = cipher.log_cap + params.log_hi_modulus
    ax = from_rns_transformed(rarot * rk.key.rax, log_modulus) >> params.log_hi_modulus
    bx = from_rns_transformed(rarot * rk.key.rbx, log_modulus) >> params.log_hi_modulus
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

    raconj = to_rns_transformed(plan, axconj, params.log_lo_modulus + params.log_hi_modulus)

    log_modulus = cipher.log_cap + params.log_hi_modulus
    ax = from_rns_transformed(raconj * ck.key.rax, log_modulus) >> params.log_hi_modulus
    bx = from_rns_transformed(raconj * ck.key.rbx, log_modulus) >> params.log_hi_modulus
    bx = bx + bxconj

    Ciphertext(
        params,
        ax,
        bx,
        cipher.log_cap,
        cipher.log_precision,
        cipher.slots)
end
