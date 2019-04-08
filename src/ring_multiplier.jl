function powmod(x::BigInt, y::Integer, m::BigInt)
    @assert y >= 0
    if y == 0
        one(T)
    elseif y == 1
        x
    else
        acc = one(BigInt)
        while y > 1
            if isodd(y)
                acc = mod(acc * x, m)
            end
            x = mod(x * x, m)
            y >>= 1
        end
        mod(x * acc, m)
    end
end



struct RingMultiplier

    pVec :: Array{UInt64, 1}
    prVec :: Array{UInt64, 1}
    pInvVec :: Array{UInt64, 1}
    scaledRootPows :: Array{Array{UInt64, 1}, 1}
    scaledRootInvPows :: Array{Array{UInt64, 1}, 1}
    scaledNInv :: Array{UInt64, 1}
    #_ntl_general_rem_one_struct* red_ss_array[nprimes];
    #mulmod_precon_t* coeffpinv_array[nprimes];

    pProd :: Array{BigInt, 1}
    pProdh :: Array{BigInt, 1}
    pHat :: Array{Array{BigInt, 1}, 1}
    pHatInvModp :: Array{Array{UInt64, 1}, 1}

    reconstruct_coeffs :: Array{BigInt, 1}
    full_modulus :: BigInt

    function RingMultiplier()

        pVec = Array{UInt64}(undef, nprimes)

        primetest = (one(UInt64) << pbnd) + one(UInt64)
        for i in 0:nprimes-1
            while true
                primetest += M
                if isprime(primetest)
                    pVec[i+1] = primetest
                    break
                end
            end
        end

        pInvVec = HEAAN_inv.(pVec)

        scaledRootPows = [Array{UInt64}(undef, N) for i in 1:nprimes]

        for i in 0:nprimes-1
            #prVec[i] = trunc(UInt64, (UInt128(1) << kbar2) ÷ pVec[i+1])
            root = findMthRootOfUnity(M, pVec[i+1])
            rootinv = invmod(root, pVec[i+1])
            NInv = invmod(N, pVec[i+1])
            #mulMod(scaledNInv[i], NInv, (1ULL << 32), pVec[i]);
            #mulMod(scaledNInv[i], scaledNInv[i], (1ULL << 32), pVec[i]);
            #scaledRootPows[i] = new uint64_t[N]();
            #scaledRootInvPows[i] = new uint64_t[N]();
            power = one(UInt64)
            powerInv = one(UInt64)
            for j in 0:N-1
                jprime = DarkIntegers.bitreverse(UInt32(j), 32) >> (32 - logN)
                rootpow = power
                scaledRootPows[i+1][jprime+1] = mulmod(rootpow, one(UInt64) << 32, pVec[i+1])
                scaledRootPows[i+1][jprime+1] = mulmod(
                    scaledRootPows[i+1][jprime+1], one(UInt64) << 32, pVec[i+1])
                #uint64_t rootpowInv = powerInv;
                #mulMod(scaledRootInvPows[i][jprime], rootpowInv, (1ULL << 32), pVec[i]);
                #mulMod(scaledRootInvPows[i][jprime], scaledRootInvPows[i][jprime], (1ULL << 32), pVec[i]);
                power = DarkIntegers.mulmod(power, root, pVec[i+1])
                powerInv = DarkIntegers.mulmod(powerInv, rootinv, pVec[i+1])
            end
        end

        prVec = Array{UInt64}(undef, nprimes)

        scaledRootInvPows = Array{Array{UInt64, 1}}(undef, nprimes)
        scaledNInv = Array{UInt64}(undef, nprimes)

        # my CRT reconstruct algorithm
        full_modulus = prod(BigInt.(pVec))
        # totient of a prime `p` is just `p-1`
        reconstruct_coeffs = [powmod(full_modulus ÷ p, p - 1, full_modulus) for p in pVec]


        pProd = Array{BigInt}(undef, nprimes)
        pProdh = Array{BigInt}(undef, nprimes)
        pHat = Array{Array{BigInt, 1}}(undef, nprimes)
        pHatInvModp = Array{Array{BigInt, 1}}(undef, nprimes)
        for i in 0:nprimes-1
            #coeffpinv_array[i] = new mulmod_precon_t[i + 1];
            pProd[i+1] = (i == 0) ? BigInt(pVec[i+1]) : (pProd[i] * pVec[i+1])
            pHat[i+1] = Array{BigInt}(undef, i+1)
            pHatInvModp[i+1] = Array{UInt64}(undef, i+1)
            for j in 0:i
                pHat[i+1][j+1] = one(BigInt)
                for k in 0:j-1
                    pHat[i+1][j+1] *= pVec[k+1]
                end
                for k in j+1:i
                    pHat[i+1][j+1] *= pVec[k+1]
                end
                pHatInvModp[i+1][j+1] = trunc(UInt64, pHat[i+1][j+1] % pVec[j+1])
                pHatInvModp[i+1][j+1] = invmod(pHatInvModp[i+1][j+1], pVec[j+1])
                #coeffpinv_array[i][j] = PrepMulModPrecon(pHatInvModp[i][j], pVec[j]);
            end
        end

        #=
        for i in 0:nprimes-1
            println("pProd: ", pProd[i+1])

            print("pHat:")
            for j in 0:i
                print(" ", pHat[i+1][j+1])
            end
            println()

            print("pHatInvModp:")
            for j in 0:i
                print(" ", pHatInvModp[i+1][j+1])
            end
            println()
        end
        exit(0)
        =#

        new(pVec, prVec, pInvVec, scaledRootPows, scaledRootInvPows, scaledNInv,
            pProd, pProdh, pHat, pHatInvModp, reconstruct_coeffs, full_modulus)
    end
end


function findMthRootOfUnity(M::Int, m::UInt64)
    tp = RRElem{UInt64, m}
    DarkIntegers.rr_value(DarkIntegers.get_root_of_one(tp, M, false))
end


function HEAAN_inv(x::UInt64)
    x^typemax(UInt64)
end


function butt(a::UInt64, b::UInt64, W::UInt64, p::UInt64, pInv::UInt64)
    U1, U0 = DarkIntegers.mulhilo(b, W)
    Q = U0 * pInv
    H, _ = DarkIntegers.mulhilo(Q, p)
    V = U1 < H ? U1 + p - H : U1 - H
    new_b = a < V ? a + p - V : a - V
    new_a = a + V
    if (new_a > p)
        new_a -= p
    end
    new_a, new_b
end


function HEAAN_NTT!(a::Array{UInt64, 1}, index::Int, rmul::RingMultiplier)
    res = copy(a)
    t = N
    logt1 = logN + 1
    p = rmul.pVec[index]
    pInv = rmul.pInvVec[index]
    for logm in 0:logN-1
        m = 1 << logm
        t >>= 1
        logt1 -= 1
        for i in 0:m-1
            j1 = i << logt1
            j2 = j1 + t - 1
            W = rmul.scaledRootPows[index][m + i + 1]
            for j in j1:j2
                n1, n2 = butt(res[j+1], res[j+t+1], W, p, pInv)
                res[j+1] = n1
                res[j+t+1] = n2
            end
        end
    end
    res
end


function multNTT(rm::RingMultiplier, a::Array{BigInt, 1}, rb::Array{UInt64, 1}, np::Int, modulus::BigInt)

    # TODO: duplicate code with CRT()

    res = Array{UInt64}(undef, np * N) # TODO technically it's always `Nnprimes`?
    for i in 1:N
        for j in 1:np
            res[(j-1)*N + i] = mod(a[i], rm.pVec[j])
        end
    end

    for j in 1:np
        r = res[(j-1)*N+1:j*N]
        m = rm.pVec[j]

        tp = RRElem{UInt64, m}
        rb_tp = tp.(rb[(j-1)*N+1:j*N])

        rr = DarkIntegers.ntt(tp.(r), inverse=false, negacyclic=true)
        rr2 = DarkIntegers.ntt(rr .* rb_tp, inverse=true, negacyclic=true)

        res[(j-1)*N+1:j*N] .= DarkIntegers.rr_value.(rr2)
    end

    # reconstruct
    res_bi = Array{BigInt}(undef, N)

    pProdnp = rm.pProd[np]
    ps = rm.pVec[1:np]
    full_modulus = prod(BigInt.(ps))
    # totient of a prime `p` is just `p-1`
    reconstruct_coeffs = [powmod(full_modulus ÷ p, p - 1, full_modulus) for p in ps]

    for i in 1:N
        r = mod(sum(BigInt.(res[i:N:end]) .* reconstruct_coeffs), full_modulus)
        res_bi[i] = mod(r > (full_modulus >> 1) ? (r - full_modulus) : r, modulus)
        #res_bi[i] = reconstruct(rm, res, i, np, modulus) # TODO: check if that is faster
    end
    res_bi
end


function reconstruct(rm::RingMultiplier, rx::Array{UInt64, 1}, n::Int, np::Int, q::BigInt)
    pHatnp = rm.pHat[np]
    pHatInvModpnp = rm.pHatInvModp[np]
    #mulmod_precon_t* coeffpinv_arraynp = coeffpinv_array[np - 1];
    pProdnp = rm.pProd[np]
    pProdhnp = pProdnp >> 1

    acc = zero(BigInt)
    for i in 0:np-1
        p = rm.pVec[i+1]
        tt = pHatInvModpnp[i+1]
        s = mulmod(rx[n + (i << logN)], tt, p)
        acc += pHatnp[i+1] * s
    end
    res = mod(acc, pProdnp)
    if res > pProdhnp
        res -= pProdnp
    end
    res = mod(res, q)
    res
end
