struct PublicKeyRNS
    rax :: RNSPolynomial
    rbx :: RNSPolynomial
end


struct EncryptionKey
    params :: Params
    key :: PublicKeyRNS

    function EncryptionKey(rng::AbstractRNG, secret_key::SecretKey)
        params = secret_key.params
        log_modulus = params.log_hi_modulus + params.log_lo_modulus
        log_plen = params.log_polynomial_length
        plen = 2^log_plen

        #ax = rand_big_int(rng, log_modulus, plen)
        ax_c = sampleUniform2(rng, log_modulus, plen)
        ax = Polynomial(BinModuloInt{BigInt, log_modulus}.(ax_c), true)

        #gg = randn(rng, plen) * params.gaussian_noise_stddev
        gg = rand_gauss(rng, plen, params.gaussian_noise_stddev)
        bx = round.(Int, gg) - secret_key * ax

        plan = rns_plan(params)
        rax = ntt_rns(to_rns(plan, ax))
        rbx = ntt_rns(to_rns(plan, bx))

        new(params, PublicKeyRNS(rax, rbx))
    end
end


struct MultiplicationKey
    params :: Params
    key :: PublicKeyRNS

    function MultiplicationKey(rng::AbstractRNG, secret_key::SecretKey)

        params = secret_key.params
        log_modulus = params.log_hi_modulus + params.log_lo_modulus
        log_plen = params.log_polynomial_length
        plen = 2^log_plen

        #ax = rand_big_int(rng, log_modulus, plen)
        ax_c = sampleUniform2(rng, log_modulus, plen)
        ax = Polynomial(BinModuloInt{BigInt, log_modulus}.(ax_c), true)

        #gg = randn(rng, plen) * params.gaussian_noise_stddev
        gg = rand_gauss(rng, plen, params.gaussian_noise_stddev)
        sxsx = square(secret_key) << params.log_hi_modulus
        bx = round.(Int, gg) - secret_key * ax + sxsx

        plan = rns_plan(params)
        rax = ntt_rns(to_rns(plan, ax))
        rbx = ntt_rns(to_rns(plan, bx))

        new(params, PublicKeyRNS(rax, rbx))
    end
end


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
        #ax = rand_big_int(rng, log_modulus, plen)
        ax_c = sampleUniform2(rng, log_modulus, plen)
        ax = Polynomial(BinModuloInt{BigInt, log_modulus}.(ax_c), true)

        #gg = randn(rng, plen) * params.gaussian_noise_stddev
        gg = rand_gauss(rng, plen, params.gaussian_noise_stddev)
        bx = round.(Int, gg) - secret_key * ax

        spow = left_rotate(as_polynomial(secret_key), shift) << params.log_hi_modulus
        bx = bx + spow

        plan = rns_plan(params)
        rax = ntt_rns(to_rns(plan, ax))
        rbx = ntt_rns(to_rns(plan, bx))

        new(params, PublicKeyRNS(rax, rbx), shift)
    end
end


struct PublicKeySet

    enc_key :: EncryptionKey
    mul_key :: MultiplicationKey
    #conj_key :: ConjugationKey
    #left_rot_keys :: Dict{Int, LeftRotationKey}

    function PublicKeySet(rng::AbstractRNG, secret_key::SecretKey)
        enc_key = EncryptionKey(rng, secret_key)
        mul_key = MultiplicationKey(rng, secret_key)
        new(enc_key, mul_key)
    end
end

