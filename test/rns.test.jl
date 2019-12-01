using HEAAN: BinModuloInt, RNSPlan, to_rns, from_rns, to_rns_signed, from_rns_signed, nprimes_for_modulus


@testgroup "RNS" begin


@testcase "full range" begin

    primes = UInt64[7, 11, 13, 17]
    plan = RNSPlan(primes)

    nprimes = 3
    p = prod(primes[1:nprimes])

    for x in 0:p-1
        x_rns = to_rns(plan, big(x), nprimes)
        x_back = from_rns(plan, x_rns)

        if x_back != x
            @test_fail "Transforming $x, got $x_back back"
            return
        end
    end
end


@testcase "binary modulus range" begin

    primes = UInt64[7, 11, 13, 17]
    plan = RNSPlan(primes)

    nprimes = 3
    p = prod(primes[1:nprimes])
    log_modulus = 9 # need 2^log_modulus <= p

    tp = BinModuloInt{BigInt, log_modulus}

    for x in 0:2^log_modulus-1
        x_tp = convert(tp, x)
        x_rns = to_rns_signed(plan, x_tp, nprimes)

        # Reference `x` transformed from all-positive to the signed range
        x_neg = is_negative(x, log_modulus) ? x - 2^log_modulus : x
        ref_rns = mod.(x_neg, primes[1:nprimes])
        if x_rns != ref_rns
            @test_fail "Transforming $x, expected $ref_rns, got $x_rns"
            return
        end

        x_back = from_rns_signed(plan, tp, x_rns)

        if x_back != x_tp
            @test_fail "Transforming $x, got $x_back back"
            return
        end
    end
end


@testcase "range selector" begin
    primes = UInt64[7, 11, 13, 17]
    plan = RNSPlan(primes)

    max_log = floor(Int, log2(prod(primes)))
    for log_modulus in 1:max_log
        nprimes = nprimes_for_modulus(plan, log_modulus)
        if prod(primes[1:nprimes]) < 2^log_modulus
            @test_fail "Got an incorrect range ($nprimes) for log_modulus=$log_modulus"
            return
        end
    end

    @test_throws DomainError nprimes_for_modulus(plan, max_log + 1)
end


end
