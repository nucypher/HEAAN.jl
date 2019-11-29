struct Params

    log_polynomial_length :: Int # in HEAAN: N
    log_lo_modulus :: Int # in HEAAN: logQ
    log_hi_modulus :: Int # in HEAAN: == logQ
    gaussian_noise_stddev :: Float64 # in HEAAN: sigma
    secret_key_length :: Int # in HEAAN: h

    function Params(; log_polynomial_length::Int=16, log_lo_modulus::Int=300)
        # TODO: this all (including log_polynomial_length) should be derived
        # from the security parameter `lambda`.

        # From Cheon et al.,
        # "We need to set the ring dimension N that satisfies the security condition
        # `N >= (λ+110) / 7.2 * (log_hi_modulus + log_lo_modulus)` ...
        # it suffices to assume that P [`log_hi_modulus`] is
        # approximately equal to qL [`log_lo_modulus`]."

        gaussian_noise_stddev = 3.2
        secret_key_length = 64
        log_hi_modulus = log_lo_modulus

        new(
            log_polynomial_length,
            log_lo_modulus,
            log_hi_modulus,
            gaussian_noise_stddev,
            secret_key_length,
            )
    end
end


function _rns_plan(params::Params)
    pbnd = 59 # TODO: move to Params?
    # TODO: why `2 * (params.log_lo_modulus + params.log_hi_modulus)`?
    # do we really use this range anywhere?
    nprimes = (
        2 + params.log_polynomial_length +
        2 * (params.log_lo_modulus + params.log_hi_modulus) + pbnd - 1) ÷ pbnd

    pVec = Array{UInt64}(undef, nprimes)
    N = 2^params.log_polynomial_length

    # TODO: precompute the primes? At least for the case of the default N
    primetest = (one(UInt64) << pbnd) + one(UInt64)
    for i in 0:nprimes-1
        while true
            primetest += N * 2
            if isprime(primetest)
                pVec[i+1] = primetest
                break
            end
        end
    end

    RNSPlan(pVec)
end


const _rns_plans = IdDict{Params, RNSPlan}()


function rns_plan(params::Params)
    if haskey(_rns_plans, params)
        _rns_plans[params]
    else
        res = _rns_plan(params)
        _rns_plans[params] = res
        res
    end
end
