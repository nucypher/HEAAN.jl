struct PublicKeyRNS
    rax :: RNSPolynomialTransformed
    rbx :: RNSPolynomialTransformed
end


"""
    EncryptionKey(rng::AbstractRNG, secret_key::SecretKey)

A public key used for encryption.
Needs a [`SecretKey`](@ref) object.
"""
struct EncryptionKey
    params :: Params
    key :: PublicKeyRNS

    function EncryptionKey(rng::AbstractRNG, secret_key::SecretKey)
        params = secret_key.params
        log_plen = params.log_polynomial_length
        plen = 2^log_plen

        # TODO: (issue #13) C++ HEAAN has `log_lo_modulus + log_hi_modulus` here
        # for the encryption key (that is, P * q_L in the notation of the paper).
        # But in both HEAAN papers only `q_L` is used as the encryption key modulus.
        log_modulus = params.log_lo_modulus + params.log_hi_modulus

        ax = rand_big_int(rng, log_modulus, plen)
        bx = discrete_gaussian(rng, params.gaussian_noise_stddev, plen) - secret_key * ax

        plan = rns_plan(params)
        rax = to_rns_transformed(plan, ax, params.log_lo_modulus)
        rbx = to_rns_transformed(plan, bx, params.log_lo_modulus)

        new(params, PublicKeyRNS(rax, rbx))
    end
end


"""
    MultiplicationKey(rng::AbstractRNG, secret_key::SecretKey)

A public key used for multiplication of two ciphertexts.
Needs a [`SecretKey`](@ref) object.
"""
struct MultiplicationKey
    params :: Params
    key :: PublicKeyRNS

    function MultiplicationKey(rng::AbstractRNG, secret_key::SecretKey)

        params = secret_key.params
        log_modulus = params.log_hi_modulus + params.log_lo_modulus
        log_plen = params.log_polynomial_length
        plen = 2^log_plen

        ax = rand_big_int(rng, log_modulus, plen)

        sxsx = broadcast_into_polynomial(
            <<, square(secret_key, params.log_lo_modulus), params.log_hi_modulus)
        bx = discrete_gaussian(rng, params.gaussian_noise_stddev, plen) - secret_key * ax + sxsx

        plan = rns_plan(params)
        rax = to_rns_transformed(plan, ax, params.log_lo_modulus)
        rbx = to_rns_transformed(plan, bx, params.log_lo_modulus)

        new(params, PublicKeyRNS(rax, rbx))
    end
end


"""
    LeftRotationKey(rng::AbstractRNG, secret_key::SecretKey, shift::Int)

A public key used for left rotation (`circshift()` with a negative argument) of a ciphertext.
Needs a [`SecretKey`](@ref) object and the absolute value of the shift
(a separate key is needed for each value of the shift).
"""
struct LeftRotationKey
    params :: Params
    key :: PublicKeyRNS
    shift :: Int

    function LeftRotationKey(rng::AbstractRNG, secret_key::SecretKey, shift::Int)
        params = secret_key.params
        log_modulus = params.log_hi_modulus + params.log_lo_modulus
        log_plen = params.log_polynomial_length
        plen = 2^log_plen

        np = cld(1 + log_modulus + log_plen + 2, 59)

        ax = rand_big_int(rng, log_modulus, plen)
        bx = discrete_gaussian(rng, params.gaussian_noise_stddev, plen) - secret_key * ax

        spow = broadcast_into_polynomial(
            <<,
            left_rotate(as_polynomial(secret_key, params.log_lo_modulus), shift),
            params.log_hi_modulus)
        bx = bx + spow

        plan = rns_plan(params)
        rax = to_rns_transformed(plan, ax, params.log_lo_modulus)
        rbx = to_rns_transformed(plan, bx, params.log_lo_modulus)

        new(params, PublicKeyRNS(rax, rbx), shift)
    end
end



"""
    LeftRotationKey(rng::AbstractRNG, secret_key::SecretKey)

A public key used for conjugation of a ciphertext.
Needs a [`SecretKey`](@ref) object.
"""
struct ConjugationKey
    params :: Params
    key :: PublicKeyRNS

    function ConjugationKey(rng::AbstractRNG, secret_key::SecretKey)
        params = secret_key.params
        log_modulus = params.log_hi_modulus + params.log_lo_modulus
        log_plen = params.log_polynomial_length
        plen = 2^log_plen

        #np = cld(1 + log_modulus + log_plen + 2, 59)

        ax = rand_big_int(rng, log_modulus, plen)
        bx = discrete_gaussian(rng, params.gaussian_noise_stddev, plen) - secret_key * ax

        sxconj = broadcast_into_polynomial(
            <<,
            conjugate(as_polynomial(secret_key, params.log_lo_modulus)),
            params.log_hi_modulus)

        bx = bx + sxconj

        plan = rns_plan(params)
        rax = to_rns_transformed(plan, ax, params.log_lo_modulus)
        rbx = to_rns_transformed(plan, bx, params.log_lo_modulus)

        new(params, PublicKeyRNS(rax, rbx))
    end
end
