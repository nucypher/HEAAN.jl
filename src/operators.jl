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


# TODO: seems to correspond to Rotate() in the paper?
# Or does it refer to `circshift` itself?
# "For an input encryption of m(Y), return an encryption of m(Y^(5^k)) in the same level"
function left_rotate(x::Polynomial, r::Integer)
    res = Polynomial(similar(x.coeffs), x.negacyclic)
    n = length(x.coeffs)
    pow = mod(5^r, 2n)
    for i in 0:n-1
        shift = mod(i * pow, 2n)
        if shift < n
            res.coeffs[shift+1] = x.coeffs[i+1]
        else
            res.coeffs[shift - n + 1] = -x.coeffs[i+1]
        end
    end
    res
end


function Base.circshift(rk::LeftRotationKey, cipher::Ciphertext, shift::Integer)
    params = cipher.params
    plan = rns_plan(params)

    # TODO: handle positive shifts too (mod N?)
    shift = -shift

    axrot = left_rotate(cipher.ax, shift)
    bxrot = left_rotate(cipher.bx, shift)

    np = cld(cipher.log_cap + params.log_lo_modulus +
        params.log_hi_modulus + params.log_polynomial_length + 2, 59)

    rarot = ntt_rns(to_rns(plan, axrot, np))

    tp_big = Polynomial{BinModuloInt{BigInt, cipher.log_cap + params.log_lo_modulus}}
    ax = from_rns(tp_big, ntt_rns(rarot * rk.key.rax, inverse=true), np) >> params.log_lo_modulus
    bx = from_rns(tp_big, ntt_rns(rarot * rk.key.rbx, inverse=true), np) >> params.log_lo_modulus
    bx = bx + bxrot

    Ciphertext(
        params,
        ax,
        bx,
        cipher.log_cap,
        cipher.log_precision,
        cipher.slots)
end
