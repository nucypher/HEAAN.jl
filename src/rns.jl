struct RNS
    pVec :: Array{UInt64, 1}
    pProd :: Array{BigInt, 1}
    pProdh :: Array{BigInt, 1}
    reconstruct_coeffs :: Array{Array{BigInt, 1}, 1}
    nprimes :: Int

    function RNS(params::Params)

        pbnd = 59 # TODO: move to Params?
        # TODO: why `2 * (params.log_lo_modulus + params.log_hi_modulus)`?
        # do we really use this range anywhere?
        nprimes = (
            2 + params.log_polynomial_length +
            2 * (params.log_lo_modulus + params.log_hi_modulus) + pbnd - 1) ÷ pbnd

        pVec = Array{UInt64}(undef, nprimes)
        N = 2^params.log_polynomial_length

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

        pVec_mp = BigInt.(pVec)
        pProd = [prod(pVec_mp[1:i]) for i in 1:nprimes]
        pProdh = pProd .>> 1

        reconstruct_coeffs = [
            [
                (pProd[i] ÷ pVec[j]) * invmod(trunc(UInt64, mod(pProd[i] ÷ pVec[j], pVec[j])), pVec[j])
                # essentially `invmod(pProd[i] ÷ pVec[j], pVec[j]) * (pProd[i] ÷ pVec[j])`,
                # but slightly faster by only calculating `invmod()` for small numbers
                for j in 1:i]
            for i in 1:nprimes]

        new(pVec, pProd, pProdh, reconstruct_coeffs, nprimes)
    end
end


function to_rns(rns::RNS, x::BigInt, np::Int, logq::Int)
    x_neg = is_negative(x, logq)
    # TODO: check that it works as intended for x == q/2 (which is negative in our definition)
    if x_neg
        x = modulus(BigInt, logq) - x
    end
    res = Array{UInt64}(undef, np)
    for j in 1:np
        p = rns.pVec[j]
        m = trunc(UInt64, mod(x, p))
        res[j] = x_neg && !iszero(m) ? p - m : m
    end
    res
end


function from_rns(rns::RNS, rx::Array{UInt64, 1}, np::Int, logq::Int)
    pProdnp = rns.pProd[np]
    pProdhnp = rns.pProdh[np]
    reconstruct_coeffs = rns.reconstruct_coeffs[np]

    acc = zero(BigInt)
    for i in 1:np
        acc += reconstruct_coeffs[i] * rx[i]
    end

    res = mod(acc, pProdnp)
    normalize(res > pProdhnp ? (res - pProdnp) : res, logq)
end


const _rns_plans = IdDict{Params, RNS}()


function rns_plan(params::Params)
    if haskey(_rns_plans, params)
        _rns_plans[params]
    else
        res = RNS(params)
        _rns_plans[params] = res
        res
    end
end
