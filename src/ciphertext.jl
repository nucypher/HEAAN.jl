struct Plaintext
    params :: Params
    polynomial :: Polynomial
    # TODO: technically, log_cap is encoded in Polynomial eltypes already
    log_cap :: Int
    log_precision :: Int
    slots :: Int
end


struct Ciphertext
    params :: Params
    ax :: Polynomial
    bx :: Polynomial
    log_cap :: Int
    log_precision :: Int
    slots :: Int

    function Ciphertext(
            params::Params, ax::Polynomial{BinModuloInt{T, Q}},
            bx::Polynomial{BinModuloInt{T, Q}},
            log_cap::Int, log_precision::Int, slots::Int) where {T, Q}
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


function rand_big_int(rng::AbstractRNG, log_modulus::Int, dims...)
    coeffs = rand(rng, zero(BigInt):((one(BigInt) << log_modulus) - one(BigInt)), dims...)
    Polynomial(BInModuloInt{BigInt, log_modulus}.(coeffs), true)
end


function encode(params::Params, vals::Array{Complex{Float64}, 1}, log_precision::Int, log_modulus::Int)
    plen = 2^params.log_polynomial_length
    slots = length(vals)
    gap = plen ÷ 2 ÷ slots
    uvals = unembed(params, vals)
    tp = BinModuloInt{BigInt, log_modulus}
    mx = zeros(tp, gap, slots, 2)
    mx[1,:,1] = float_to_integer.(tp, real.(uvals), log_precision)
    mx[1,:,2] = float_to_integer.(tp, imag.(uvals), log_precision)
    poly = Polynomial(mx[:], true)
    Plaintext(params, poly, log_modulus, log_precision, slots)
end


function bit(x::BigInt, i::Int)
    x & (one(BigInt) << i) != 0
end

function sample_ZO(rng, len::Int, log_modulus::Int)
    res = Array{BigInt}(undef, len)
    tmp = myRandomBits_ZZ(rng, len*2)
    for i in 0:len-1
        res[i+1] = (!bit(tmp, 2 * i)) ? zero(BigInt) : (!bit(tmp, 2 * i + 1)) ? one(BigInt) : -one(BigInt)
    end
    Polynomial(BinModuloInt{BigInt, log_modulus}.(normalize.(res, log_modulus)), true)
end

#=
function sample_ZO(rng::AbstractRNG, len::Int, log_modulus::Int)
    modulus = one(BigInt) << log_modulus
    minus_one = modulus - one(BigInt)

    bits1 = rand(rng, Bool, len)
    bits2 = rand(rng, Bool, len)
    res = [
        b1 ? zero(BigInt) : (b2 ? one(BigInt) : minus_one)
        for (b1, b2) in zip(bits1, bits2)]

    CappedPolynomial(Polynomial(res, true), log_modulus)
end
=#

function encrypt(rng::AbstractRNG, key::EncryptionKey, plain::Plaintext)

    params = plain.params
    plen = 2^params.log_polynomial_length

    log_modulus = plain.log_cap
    modulus = one(BigInt) << log_modulus

    vx = sample_ZO(rng, 2^params.log_polynomial_length, log_modulus)

    # TODO: `log_modulus` should be enough here?
    np = cld(1 + params.log_lo_modulus + params.log_hi_modulus + params.log_polynomial_length + 2, 59)

    # TODO: move round(...) * ... into a function of Params object
    # gg = randn(rng, plen) * params.gaussian_noise_stddev
    gg = rand_gauss(rng, plen, params.gaussian_noise_stddev)
    ax = (
        round.(Int, gg) +
        mult(vx, key.key.rax, np))

    # gg = randn(rng, plen) * params.gaussian_noise_stddev
    gg = rand_gauss(rng, plen, params.gaussian_noise_stddev)
    bx = (
        round.(Int, gg) +
        plain.polynomial +
        mult(vx, key.key.rbx, np))

    ax = ax >> params.log_lo_modulus
    bx = bx >> params.log_lo_modulus

    Ciphertext(params, ax, bx,
        plain.log_cap - params.log_lo_modulus,
        plain.log_precision - params.log_lo_modulus, plain.slots)
end


function encrypt(
        rng::AbstractRNG, key::EncryptionKey, vals::Array{Complex{Float64}, 1},
        log_precision::Int, log_cap::Int)
    params = key.params
    plain = encode(
        params, vals,
        log_precision + params.log_lo_modulus, log_cap + params.log_lo_modulus)
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


function decrypt(secret_key::SecretKey, cipher::Ciphertext)
    plain = decrypt_to_plaintext(secret_key, cipher)
    decode(plain)
end
