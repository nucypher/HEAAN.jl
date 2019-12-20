"""
    Params(; log_polynomial_length::Int=16, log_lo_modulus::Int=300)

HEAAN scheme parameters.

`log_polynomial_length` determines the length of polynomials
used for keys and ciphertexts; larger length is more secure, but also slower.

`log_lo_modulus` determines the maximum computation resource that a ciphertext can have.
Larger values are more secure and allow for longer computation sequences without bootstrapping;
but, again, are slower and make ciphertexts take more place in memory.
"""
struct Params

    log_polynomial_length :: Int # in HEAAN: N
    log_lo_modulus :: Int # in HEAAN: logQ, q_L in the paper
    log_hi_modulus :: Int # in HEAAN: == logQ, P in the paper
    gaussian_noise_stddev :: Float64 # in HEAAN: sigma
    secret_key_length :: Int # in HEAAN: h

    function Params(; log_polynomial_length::Int=16, log_lo_modulus::Int=300)
        # TODO: (issue #14) this all (including log_polynomial_length) should be derived
        # from the security parameter `lambda`.

        # From Cheon et al.,
        # "We need to set the ring dimension N that satisfies the security condition
        # `N >= (λ+110) / 7.2 * (log_hi_modulus + log_lo_modulus)` ...
        # it suffices to assume that P [`log_hi_modulus`] is
        # approximately equal to qL [`log_lo_modulus`]."

        gaussian_noise_stddev = 3.2
        secret_key_length = 64
        log_hi_modulus = log_lo_modulus + 50 # To detect misplaced lo/hi modulus

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
    # Maximum range for RNS numbers
    # (lo_modulus + hi_modulus) is the maximum key size, which will be multiplied by
    # lo_modulus-sized polynomial, which will require an additional log_polynomial_length range.
    # There may be +1 or +2 necessary, but in any case it will be caught by assertions in rns.jl.
    log_max_range = (
        params.log_polynomial_length +
        params.log_lo_modulus + params.log_hi_modulus +
        params.log_lo_modulus)
    max_range = one(BigInt) << log_max_range

    primes_prod = one(BigInt)
    primes = UInt64[]
    N = 2^params.log_polynomial_length

    # In practice, the primes can be precomputed,
    # but for our parameters it takes very little time (and the plan is cached anyway),
    # so for simplicity we don't bother.
    #
    # We're only using primes with 62 bits or less, to make addmod() a little faster
    # (since there's no UInt64 overflow). Otherwise we could start from 2^64-1.
    primetest = (one(UInt64) << 63) + one(UInt64)
    while true
        primetest -= N * 2 # need (prime - 1) to be a multiple of 2N for NTT to work
        if isprime(primetest)
            push!(primes, primetest)
            primes_prod *= big(primetest)
            if primes_prod > max_range
                break
            end
        end
    end

    RNSPlan(primes)
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


const _embedding_plans = IdDict{Params, EmbeddingPlan}()


function embedding_plan(params::Params)
    if haskey(_embedding_plans, params)
        _embedding_plans[params]
    else
        res = EmbeddingPlan(2^params.log_polynomial_length ÷ 2)
        _embedding_plans[params] = res
        res
    end
end


function embed(params::Params, vals::Array{Complex{Float64}, 1})
    embed(embedding_plan(params), vals)
end


function unembed(params::Params, vals::Array{Complex{Float64}, 1})
    unembed(embedding_plan(params), vals)
end
