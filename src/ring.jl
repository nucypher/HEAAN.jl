struct BootContext
end


struct Ring

    qpows :: Array{BigInt, 1}
    rotGroup :: Array{Int, 1}
    ksiPows :: Array{Complex{Float64}, 1}
    bootContextMap :: Dict{Int, BootContext}
    multiplier :: RingMultiplier

    function Ring()

        qpows = Array{BigInt}(undef, logQQ + 1)
        qpows[1] = one(BigInt)
        for i in 2:logQQ+1
            qpows[i] = qpows[i - 1] << 1
        end

        rotGroup = Array{Int}(undef, Nh)
        fivePows = 1
        for i in 1:Nh
            rotGroup[i] = fivePows
            fivePows *= 5
            fivePows = mod(fivePows, M)
        end

        ksiPows = Array{Complex{Float64}}(undef, M + 1)
        m_pi = Float64(pi)
        for j in 1:M
            angle = 2.0 * m_pi * (j-1) / M
            ksiPows[j] = mycos(angle) + im * mysin(angle)
        end
        ksiPows[M + 1] = ksiPows[1]

        new(qpows, rotGroup, ksiPows, Dict{Int, BootContext}(), RingMultiplier())
    end

end


function sampleHWT(rng::MyRNG, ring::Ring)
    res = zeros(BigInt, N)
    idx = 0
    tmp = myRandomBits_ZZ(rng, h)
    #println("Sampled: $tmp")
    while idx < h
        i = myRandomBits_long(rng, logN)
        if res[i+1] == 0
            res[i+1] = (tmp & (one(BigInt) << idx) == 0) ? one(BigInt) : (-one(BigInt))
            #println("Filling position $i with $(res[i+1])")
            idx += 1
        end
    end
    res
end


function sampleUniform2(rng::MyRNG, bits::Int)
    [myRandomBits_ZZ(rng, bits) for i in 1:N]
end


function zero_high_bits(x::MPNumber{N, T}, l::Int) where {N, T}
    if l == 0
        x
    else
        Base.setindex(x, (x[N] << l) >> l, N)
    end
end


function mult(ring::Ring, a::Array{BigInt, 1}, b::Array{BigInt, 1}, np::Int, q::BigInt)

    # TODO: technically `a` and `b` can have negative coefficients

    @assert (one(BigInt) << trailing_zeros(q)) == q
    logq = trailing_zeros(q)
    n = cld(logq, 64)
    tp = MPNumber{n, UInt64}

    aa = Polynomial{tp}(a, true)
    bb = Polynomial{tp}(b, true)
    res = (aa * bb).coeffs

    res .= zero_high_bits.(res, n * 64 - logq)

    convert.(BigInt, res)

    #mod.((Polynomial{BigInt}(a, true) * Polynomial{BigInt}(b, true)).coeffs, q)
end


function subFromGaussAndEqual(rng::MyRNG, ring::Ring, x::Array{BigInt, 1}, q::BigInt)

    bignum = Int(0xfffffff)

    res = similar(x)

    for i in 0:2:N-1
        # TODO: essentially Box-Muller sampling
        r1 = (1 + myrand_long(rng, bignum)) / (bignum + 1)
        r2 = (1 + myrand_long(rng, bignum)) / (bignum + 1)
        theta = 2pi * r1
        rr = sqrt(-2.0 * log(r2)) * sigma

        res[i+1] = AddMod(-x[i+1], floor(BigInt, rr * cos(theta) + 0.5), q)
        res[i+2] = AddMod(-x[i+2], floor(BigInt, rr * sin(theta) + 0.5), q)
    end

    res
end


function addGaussAndEqual(rng, ring::Ring, x::Array{BigInt, 1}, q::BigInt)

    bignum = Int(0xfffffff)
    res = similar(x)

    for i in 0:2:N-1
        # TODO: essentially Box-Muller sampling
        r1 = (1 + myrand_long(rng, bignum)) / (bignum + 1)
        r2 = (1 + myrand_long(rng, bignum)) / (bignum + 1)
        theta = 2pi * r1
        rr = sqrt(-2.0 * log(r2)) * sigma

        res[i+1] = AddMod(x[i+1], floor(BigInt, rr * cos(theta) + 0.5), q)
        res[i+2] = AddMod(x[i+2], floor(BigInt, rr * sin(theta) + 0.5), q)
    end
    res
end


function CRT(ring::Ring, x::Array{BigInt, 1}, np::Int)

    res = Array{UInt64}(undef, np * N) # TODO technically it's always `Nnprimes`?
    for i in 1:N
        for j in 1:np
            res[(j-1)*N + i] = mod(x[i], ring.multiplier.pVec[j])
        end
    end

    for j in 1:np
        r = res[(j-1)*N+1:j*N]
        m = ring.multiplier.pVec[j]
        tp = RRElem{UInt64, m}
        rr = DarkIntegers.ntt(tp.(r), negacyclic=true)

        res[(j-1)*N+1:j*N] .= DarkIntegers.rr_value.(rr)
    end

    res
end


function leftShiftAndEqual!(p::Array{BigInt, 1}, bits::Int, modulus::BigInt)
    for i in 0:N-1
        p[i+1] <<= bits
        p[i+1] %= modulus
    end
end


function addAndEqual!(p1::Array{BigInt, 1}, p2::Array{BigInt, 1}, modulus::BigInt)
    for i in 0:N-1
        p1[i+1] = AddMod(p1[i+1], p2[i+1], modulus)
    end
end


function sampleZO(rng, ring::Ring)
    res = Array{BigInt}(undef, N)
    tmp = myRandomBits_ZZ(rng, M)
    for i in 0:N-1
        res[i+1] = (!bit(tmp, 2 * i)) ? zero(BigInt) : (!bit(tmp, 2 * i + 1)) ? one(BigInt) : -one(BigInt)
    end
    res
end


function rightShiftAndEqual!(p::Array{BigInt, 1}, bits::Int)
    for i in 0:N-1
        p[i+1] >>= bits
    end
end


function arrayBitReverse!(vals::Array{Complex{Float64}, 1}, n::Int)
    j = 0
    for i in 1:n-1
        bit = n >> 1
        while j >= bit
            j -= bit
            bit >>= 1
        end
        j += bit
        if i < j
            t = vals[i+1]
            vals[i+1] = vals[j+1]
            vals[j+1] = t
        end
    end
end


function EMBInvLazy(ring::Ring, vals::Array{Complex{Float64}, 1}, n::Int)
    # TODO: assuming n is a power of 2
    log_n = trailing_zeros(n)
    new_vals = copy(vals)

    for log_len in log_n:-1:0
        len = 2^log_len
        for i in 0:len:n-1
            lenh = len >> 1
            lenq = len << 2
            gap = M ÷ lenq
            for j in 0:lenh-1
                idx = (lenq - (ring.rotGroup[j+1] % lenq)) * gap
                u = new_vals[i + j + 1] + new_vals[i + j + lenh + 1]
                v = new_vals[i + j + 1] - new_vals[i + j + lenh + 1]
                v *= ring.ksiPows[idx + 1]
                new_vals[i + j + 1] = u
                new_vals[i + j + lenh + 1] = v
            end
        end
    end

    arrayBitReverse!(new_vals, n)
    new_vals
end


function EMBInv(ring::Ring, vals::Array{Complex{Float64}, 1}, n::Int)
    vals = EMBInvLazy(ring, vals, n)
    for i in 0:n-1
        vals[i+1] /= n
    end
    vals
end


function EMB(ring::Ring, vals::Array{Complex{Float64}, 1}, n::Int)
    new_vals = copy(vals)
    arrayBitReverse!(new_vals, n)

    # TODO: assuming n is a power of 2
    log_n = trailing_zeros(n)

    for log_len in 1:log_n
        len = 1 << log_len

        for i in 0:len:n-1
            lenh = len >> 1
            lenq = len << 2
            gap = M ÷ lenq
            for j in 0:lenh-1
                idx = (ring.rotGroup[j+1] % lenq) * gap
                u = new_vals[i + j + 1]
                v = new_vals[i + j + lenh + 1]
                v *= ring.ksiPows[idx + 1]
                new_vals[i + j + 1] = u + v
                new_vals[i + j + lenh + 1] = u - v
            end
        end
    end

    new_vals
end


function encode(ring::Ring, mx::Array{BigInt, 1}, vals::Array{Complex{Float64}, 1}, slots::Int, logp::Int)
    gap = Nh ÷ slots
    uvals = EMBInv(ring, vals, slots)
    for i in 0:slots-1
        jdx = Nh + i * gap
        idx = i * gap
        mx[idx+1] = scaleUpToZZ(real(uvals[i+1]), logp)
        mx[jdx+1] = scaleUpToZZ(imag(uvals[i+1]), logp)
    end
end


function decode(ring::Ring, mx::Array{BigInt, 1}, slots::Int, logp::Int, logq::Int)

    vals = Array{Complex{Float64}}(undef, slots)

    q = ring.qpows[logq+1]
    gap = Nh ÷ slots

    for i in 0:slots-1
        idx = i * gap
        tmp = rem(mx[idx+1], q)
        if NumBits(tmp) == logq
            tmp -= q
        end
        vr = scaleDownToReal(tmp, logp)

        tmp = rem(mx[idx + Nh + 1], q)
        if NumBits(tmp) == logq
            tmp -= q
        end
        vi = scaleDownToReal(tmp, logp)

        vals[i+1] = vr + im * vi
    end
    EMB(ring, vals, slots)
end
