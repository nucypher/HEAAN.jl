#=
struct PublicKeySet

    # contain Encryption, Multiplication and Conjugation keys, if generated
    keyMap :: Dict{Int, Key}

    # contain left rotation keys, if generated
    #leftRotKeyMap :: Dict{Int, Key}

    # contain Encryption, Multiplication and Conjugation keys, if generated
    #serKeyMap :: Dict{Int, String}

    # contain left rotation keys, if generated
    #serLeftRotKeyMap :: Dict{Int, String}


    function PublicKeySet(rng::AbstractRNG, params::Params, secret_key::Secretkey)

    end

end
=#


struct PublicKey{T <: CappedPolynomial}
    ax :: T
    bx :: T
end


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
        ax = CappedPolynomial(Polynomial(ax_c, true), log_modulus)

        #gg = randn(rng, plen) * params.gaussian_noise_stddev
        gg = rand_gauss(rng, plen, params.gaussian_noise_stddev)
        bx = round.(Int, gg) - secret_key * ax

        plan = rns_plan(params)
        # TODO: `plan.nprimes` should be the default?
        rax = ntt_rns(to_rns(plan, ax, plan.nprimes))
        rbx = ntt_rns(to_rns(plan, bx, plan.nprimes))

        new(params, PublicKeyRNS(rax, rbx))
    end
end