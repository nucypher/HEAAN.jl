struct Plaintext
    params :: Params
    polynomial :: Polynomial
    # TODO: (issue #11) technically, log_cap is encoded in Polynomial eltypes already
    log_cap :: Int
    log_precision :: Int
    slots :: Int
end


"""
A ciphertext object created by [`encrypt`](@ref).
"""
struct Ciphertext
    params :: Params
    ax :: Polynomial
    bx :: Polynomial

    "The limit for `log_precision`, so essentially the available \"computation resource.\""
    log_cap :: Int

    """
    The amount of ciphertext occupied by the actual value. Cannot get larger than `log_cap`.
    Generally increased during arithmetic operations.
    """
    log_precision :: Int
    slots :: Int

    function Ciphertext(
            params::Params, ax::Polynomial{BinModuloInt{T, Q}},
            bx::Polynomial{BinModuloInt{T, Q}},
            log_cap::Int, log_precision::Int, slots::Int) where {T, Q}
        @assert log_cap <= params.log_lo_modulus
        @assert log_cap == Q
        @assert log_precision > 0
        @assert log_precision <= Q
        new(params, ax, bx, log_cap, log_precision, slots)
    end
end


function compatible(c1::Ciphertext, c2::Ciphertext; different_precision::Bool=false)
    (c1.params == c2.params &&
    c1.log_cap == c2.log_cap &&
    (different_precision || c1.log_precision == c2.log_precision) &&
    c1.slots == c2.slots)
end


function encode(
        params::Params, vals::Array{Complex{Float64}, 1}, log_precision::Int, log_modulus::Int,
        minimize_range::Bool=false)
    plen = 2^params.log_polynomial_length
    slots = length(vals)
    gap = plen ÷ 2 ÷ slots
    uvals = unembed(params, vals)
    tp = BinModuloInt{BigInt, log_modulus}
    mx = zeros(tp, gap, slots, 2)
    mx[1,:,1] = float_to_integer.(tp, real.(uvals), log_precision)
    mx[1,:,2] = float_to_integer.(tp, imag.(uvals), log_precision)
    poly = Polynomial(mx[:], negacyclic_modulus)

    if minimize_range
        nb = maximum(num_bits.(mx))
        # +1 to make sure all numbers from `-m` to `m` fit into the type
        # (where `m` is the maximum number that has `num_bits == nb`, that is `2^(nb+1)-1`).
        new_log_modulus = nb + 1
        poly = mod_down_to(poly, new_log_modulus)
    end

    Plaintext(params, poly, log_modulus, log_precision, slots)
end


function encrypt(rng::AbstractRNG, key::EncryptionKey, plain::Plaintext)

    params = plain.params
    plen = 2^params.log_polynomial_length

    log_modulus = plain.log_cap

    vx = sample_ZO(rng, 2^params.log_polynomial_length)

    plan = rns_plan(params)
    rvx = to_rns_transformed(plan, vx, params.log_lo_modulus + params.log_hi_modulus)

    ax = (
        discrete_gaussian(rng, params.gaussian_noise_stddev, plen) +
        from_rns_transformed(rvx * key.key.rax, log_modulus))

    bx = (
        discrete_gaussian(rng, params.gaussian_noise_stddev, plen) +
        plain.polynomial +
        from_rns_transformed(rvx * key.key.rbx, log_modulus))

    ax = broadcast_into_polynomial(right_shift_rounded, ax, params.log_hi_modulus)
    bx = broadcast_into_polynomial(right_shift_rounded, bx, params.log_hi_modulus)

    Ciphertext(params, ax, bx,
        plain.log_cap - params.log_hi_modulus,
        plain.log_precision - params.log_hi_modulus, plain.slots)
end


"""
    encrypt(
        rng::AbstractRNG, key::EncryptionKey, vals::Array{Complex{Float64}, 1},
        log_precision::Int, log_cap::Int)

Creates a [`Ciphertext`](@ref) object out of a complex-valued array `vals`.

`log_precision` determines the (absolute) precision used for encoding
(cannot be larger than `log_cap`).

`log_cap` is the available computation resource
(cannot be larger than [`Params.log_lo_modulus`](@ref Params); see explanation thereof).
"""
function encrypt(
        rng::AbstractRNG, key::EncryptionKey, vals::Array{Complex{Float64}, 1},
        log_precision::Int, log_cap::Int)
    params = key.params
    plain = encode(
        params, vals,
        log_precision + params.log_hi_modulus, log_cap + params.log_hi_modulus)
    encrypt(rng, key, plain)
end


function decrypt_to_plaintext(secret_key::SecretKey, cipher::Ciphertext)
    params = secret_key.params
    polynomial = secret_key * cipher.ax + cipher.bx
    Plaintext(params, polynomial, cipher.log_cap, cipher.log_precision, cipher.slots)
end


function decode(plain::Plaintext)

    slots = plain.slots
    vals = Array{Complex{Float64}}(undef, slots)

    log_modulus = plain.log_cap
    mx = plain.polynomial.coeffs

    Nh = 2^plain.params.log_polynomial_length ÷ 2
    gap = Nh ÷ slots

    for i in 0:slots-1
        idx = i * gap
        vr = integer_to_float(Float64, mx[idx+1], plain.log_precision)
        vi = integer_to_float(Float64, mx[idx + Nh + 1], plain.log_precision)
        vals[i+1] = vr + im * vi
    end
    embed(plain.params, vals)
end


"""
    decrypt(secret_key::SecretKey, cipher::Ciphertext)

Decrypts a [`Ciphertext`](@ref) object and returns the resulting complex-valued array.
"""
function decrypt(secret_key::SecretKey, cipher::Ciphertext)
    plain = decrypt_to_plaintext(secret_key, cipher)
    decode(plain)
end


function mul_by_rns(cipher::Ciphertext, p::RNSPolynomialTransformed, log_precision::Int)
    Ciphertext(
        cipher.params,
        mul_by_rns(cipher.ax, p),
        mul_by_rns(cipher.bx, p),
        cipher.log_cap,
        cipher.log_precision + log_precision,
        cipher.slots)
end
