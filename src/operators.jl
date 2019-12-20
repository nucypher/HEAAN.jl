"""
    add(cipher1::Ciphertext, cipher2::Ciphertext)

Elementwise addition of encrypted vectors.
Both ciphertexts must have the same `log_cap` and `log_precision`.
"""
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


"""
    sub(cipher1::Ciphertext, cipher2::Ciphertext)

Elementwise subtraction of encrypted vectors.
Both ciphertexts must have the same `log_cap` and `log_precision`.
"""
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


"""
    mul(key::MultiplicationKey, cipher1::Ciphertext, cipher2::Ciphertext)

Elementwise multiplication of encrypted vectors.
Needs a [`MultiplicationKey`](@ref).
Both ciphertexts must have the same `log_cap` (but can have different `log_precision`);
`log_precision` of the result is the sum of `log_precision` of both ciphertexts.
"""
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
    ax = broadcast_into_polynomial(
        right_shift_rounded,
        from_rns_transformed(raa * key.rax, log_modulus), params.log_hi_modulus)
    bx = broadcast_into_polynomial(
        right_shift_rounded,
        from_rns_transformed(raa * key.rbx, log_modulus), params.log_hi_modulus)

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


"""
    square(mk::MultiplicationKey, cipher::Ciphertext)

Elementwise square of an encrypted vector.
Needs a [`MultiplicationKey`](@ref) object.
`log_precision` of the result is doubled.
"""
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
    ax = broadcast_into_polynomial(
        right_shift_rounded,
        from_rns_transformed(raa * key.rax, log_modulus), params.log_hi_modulus)
    bx = broadcast_into_polynomial(
        right_shift_rounded,
        from_rns_transformed(raa * key.rbx, log_modulus), params.log_hi_modulus)

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


"""
    negate(cipher::Ciphertext)

Elementwise negation of an encrypted vector.
"""
function negate(cipher::Ciphertext)
    Ciphertext(
        cipher.params,
        -cipher.ax,
        -cipher.bx,
        cipher.log_cap,
        cipher.log_precision,
        cipher.slots)
end


"""
    mod_down_by(cipher::Ciphertext, dlog_cap::Int)

Decreases `log_cap` of a [`Ciphertext`](@ref) object by `dlog_cap`.
"""
function mod_down_by(cipher::Ciphertext, dlog_cap::Int)
    Ciphertext(
        cipher.params,
        mod_down_by(cipher.ax, dlog_cap),
        mod_down_by(cipher.bx, dlog_cap),
        cipher.log_cap - dlog_cap,
        cipher.log_precision,
        cipher.slots)
end


"""
    mod_down_to(cipher::Ciphertext, log_cap::Int)

Decreases `log_cap` of a [`Ciphertext`](@ref) object to the provided value.
"""
function mod_down_to(cipher::Ciphertext, log_cap::Int)
    Ciphertext(
        cipher.params,
        mod_down_to(cipher.ax, log_cap),
        mod_down_to(cipher.bx, log_cap),
        log_cap,
        cipher.log_precision,
        cipher.slots)
end


"""
    rescale_by(cipher::Ciphertext, dlog_cap::Int)

Decreases `log_cap` and `log_precision` of a [`Ciphertext`](@ref) object by `dlog_cap`.
"""
function rescale_by(cipher::Ciphertext, dlog_cap::Int)
    Ciphertext(
        cipher.params,
        broadcast_into_polynomial(right_shift_rounded, cipher.ax, dlog_cap),
        broadcast_into_polynomial(right_shift_rounded, cipher.bx, dlog_cap),
        cipher.log_cap - dlog_cap,
        cipher.log_precision - dlog_cap,
        cipher.slots)
end


"""
    div_by_po2(cipher::Ciphertext, bits::Int)

Elementwise division of an encrypted vector by `2^bits`.
"""
function div_by_po2(cipher::Ciphertext, bits::Int)
    Ciphertext(
        cipher.params,
        broadcast_into_polynomial(right_shift_rounded, cipher.ax, bits),
        broadcast_into_polynomial(right_shift_rounded, cipher.bx, bits),
        cipher.log_cap - bits,
        cipher.log_precision,
        cipher.slots)
end


"""
    add_const(cipher::Ciphertext, cnst::Float64)

Adds a floating-point number `cnst` to each element of an encrypted vector.
Ciphertext's `log_precision` is used to encode `cnst`.
"""
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


"""
    add_const(cipher::Ciphertext, cnst::Float64)

Adds a complex number `cnst` to each element of an encrypted vector.
Ciphertext's `log_precision` is used to encode `cnst`.
"""
function add_const(cipher::Ciphertext, cnst::Complex{Float64})
    negate(imul(add_const(imul(add_const(cipher, real(cnst))), -imag(cnst))))
end


"""
    mul_by_const(cipher::Ciphertext, cnst::Float64, log_precision::Int)

Multiplies each element of an encrypted vector by a floating-point number `cnst`.
The provided `log_precision` is used to encode `cnst`.
"""
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


"""
    mul_by_const(cipher::Ciphertext, cnst::Complex{Float64}, log_precision::Int)

Multiplies each element of an encrypted vector by a complex number `cnst`.
The provided `log_precision` is used to encode `cnst`.
"""
function mul_by_const(cipher::Ciphertext, cnst::Complex{Float64}, log_precision::Int)
    add(
        mul_by_const(cipher, real(cnst), log_precision),
        imul(mul_by_const(cipher, imag(cnst), log_precision)))
end


"""
    imul(cipher::Ciphertext)

Multiplies each element of an encrypted vector by the imaginary unit.
"""
function imul(cipher::Ciphertext)
    params = cipher.params
    shift = 1 << (params.log_polynomial_length - 1)
    Ciphertext(
        params,
        mul_by_monomial(cipher.ax, shift),
        mul_by_monomial(cipher.bx, shift),
        cipher.log_cap,
        cipher.log_precision,
        cipher.slots)
end


"""
    mul_by_const_vec(cipher::Ciphertext, v::Array{Complex{Float64}, 1}, log_precision::Int)

Multiplies each element of an encrypted vector by the corresponding element of `v`.
`log_precision` is used to encode elements of `v`.
"""
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


"""
    circshift(rk::LeftRotationKey, cipher::Ciphertext, shift::Integer)

Rotates an encrypted vector (currently only to the left, that is `shift` must be negative).
The [`LeftRotationKey`](@ref) must correspond to `shift`.
"""
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
    ax = broadcast_into_polynomial(
        right_shift_rounded,
        from_rns_transformed(rarot * rk.key.rax, log_modulus), params.log_hi_modulus)
    bx = broadcast_into_polynomial(
        right_shift_rounded,
        from_rns_transformed(rarot * rk.key.rbx, log_modulus), params.log_hi_modulus)
    bx = bx + bxrot

    Ciphertext(
        params,
        ax,
        bx,
        cipher.log_cap,
        cipher.log_precision,
        cipher.slots)
end


"""
    conj(ck::ConjugationKey, cipher::Ciphertext)

Conjugates an encrypted vector.
Needs a [`ConjugationKey`](@ref) object.
"""
function Base.conj(ck::ConjugationKey, cipher::Ciphertext)
    params = cipher.params
    plan = rns_plan(params)

    axconj = conjugate(cipher.ax)
    bxconj = conjugate(cipher.bx)

    raconj = to_rns_transformed(plan, axconj, params.log_lo_modulus + params.log_hi_modulus)

    log_modulus = cipher.log_cap + params.log_hi_modulus
    ax = broadcast_into_polynomial(
        right_shift_rounded,
        from_rns_transformed(raconj * ck.key.rax, log_modulus), params.log_hi_modulus)
    bx = broadcast_into_polynomial(
        right_shift_rounded,
        from_rns_transformed(raconj * ck.key.rbx, log_modulus), params.log_hi_modulus)
    bx = bx + bxconj

    Ciphertext(
        params,
        ax,
        bx,
        cipher.log_cap,
        cipher.log_precision,
        cipher.slots)
end
