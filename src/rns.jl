struct RNSPlan
    primes :: Array{UInt64, 1}
    primes_prod :: Array{BigInt, 1}
    primes_prod_half :: Array{BigInt, 1}
    reconstruct_coeffs :: Array{Array{BigInt, 1}, 1}
    max_bin_moduli :: Array{Int, 1}

    function RNSPlan(primes::Array{UInt64, 1})
        nprimes = length(primes)
        primes_mp = big.(primes)
        primes_prod = [prod(primes_mp[1:i]) for i in 1:nprimes]
        primes_prod_half = primes_prod .>> 1

        reconstruct_coeffs = [
            [
                (primes_prod[i] ÷ primes[j]) *
                    invmod(trunc(UInt64, mod(primes_prod[i] ÷ primes[j], primes[j])), primes[j])
                # essentially `invmod(primes_prod[i] ÷ primes[j], primes[j]) *
                #    (primes_prod[i] ÷ primes[j])`,
                # but slightly faster by only calculating `invmod()` for small numbers
                for j in 1:i]
            for i in 1:nprimes]

        max_bin_moduli = [num_bits(pp) - 1 for pp in primes_prod]

        new(primes, primes_prod, primes_prod_half, reconstruct_coeffs, max_bin_moduli)
    end
end


function nprimes_for_modulus(plan::RNSPlan, log_modulus::Int)
    idx = findfirst(x -> x >= log_modulus, plan.max_bin_moduli)
    if idx === nothing
        throw(DomainError(
            log_modulus, "Not enough primes in the plan to support this log2(modulus)"))
    end
    idx
end


function to_rns(plan::RNSPlan, x::BigInt, np::Int)
    trunc.(UInt64, mod.(x, (@view plan.primes[1:np])))
end


function negate(x::UInt64, p::UInt64)
    iszero(x) ? x : p - x
end


function to_rns(plan::RNSPlan, x::BigInt, np::Int, log_modulus::Int)
    x_neg = is_negative(x, log_modulus)
    if x_neg
        x = modulus(BigInt, log_modulus) - x
    end
    res = to_rns(plan, x, np)
    if x_neg
        negate.(res, (@view plan.primes[1:np]))
    else
        res
    end
end


function from_rns(plan::RNSPlan, rx::Array{UInt64, 1})
    np = length(rx)
    reconstruct_coeffs = plan.reconstruct_coeffs[np]

    acc = zero(BigInt)
    for i in 1:np
        acc += reconstruct_coeffs[i] * rx[i]
    end
    mod(acc, plan.primes_prod[np])
end


function from_rns(plan::RNSPlan, rx::Array{UInt64, 1}, log_modulus::Int)
    np = length(rx)
    primes_prod_half = plan.primes_prod_half[np]
    primes_prod = plan.primes_prod[np]
    res = from_rns(plan, rx)
    normalize(res > primes_prod_half ? (res - primes_prod) : res, log_modulus)
end
