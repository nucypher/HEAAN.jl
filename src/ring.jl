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
    q = one(BigInt) << logQQ
    while idx < h
        i = myRandomBits_long(rng, logN)
        if res[i+1] == 0
            res[i+1] = (tmp & (one(BigInt) << idx) == 0) ? one(BigInt) : (q - one(BigInt))
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


function zero_high_bits(x::BigInt, l::Int) where {N, T}
    @assert !signbit(x)
    x & ((one(BigInt) << l) - one(BigInt))
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

    # x is in range [0, q)

    bignum = Int(0xfffffff)

    res = similar(x)

    for i in 0:2:N-1
        # TODO: essentially Box-Muller sampling
        r1 = (1 + myrand_long(rng, bignum)) / (bignum + 1)
        r2 = (1 + myrand_long(rng, bignum)) / (bignum + 1)
        theta = 2pi * r1
        rr = sqrt(-2.0 * log(r2)) * sigma

        y1 = floor(BigInt, rr * mycos(theta) + 0.5)
        y2 = floor(BigInt, rr * mysin(theta) + 0.5)

        if signbit(y1)
            y1 = q + y1
        end
        if signbit(y2)
            y2 = q + y2
        end

        x1 = NegMod(x[i+1], q)
        x2 = NegMod(x[i+2], q)

        res[i+1] = AddMod(x1, y1, q)
        res[i+2] = AddMod(x2, y2, q)
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

        res[i+1] = AddMod(x[i+1], floor(BigInt, rr * mycos(theta) + 0.5), q)
        res[i+2] = AddMod(x[i+2], floor(BigInt, rr * mysin(theta) + 0.5), q)
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
    logq = trailing_zeros(modulus)
    for i in 0:N-1
        p[i+1] <<= bits
        p[i+1] = zero_high_bits(p[i+1], logq)
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
    tmp = one(BigInt) << (bits - 1)
    for i in 0:N-1
        p[i+1] += tmp
        if p[i+1] >= 0
            p[i+1] = p[i+1] >> bits
        else
            # TODO: mimicking NTL (and C++?) behavior here.
            # perhaps, can be removed if all the numbers are made positive.
            p[i+1] = -((-p[i+1]) >> bits)
        end
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


function encode(ring::Ring, mx::Array{BigInt, 1}, vals::Array{Complex{Float64}, 1}, slots::Int, logp::Int, log_full::Int)
    gap = Nh ÷ slots
    uvals = EMBInv(ring, vals, slots)
    for i in 0:slots-1
        jdx = Nh + i * gap
        idx = i * gap
        mx[idx+1] = float_to_integer(real(uvals[i+1]), logp, log_full)
        mx[jdx+1] = float_to_integer(imag(uvals[i+1]), logp, log_full)
    end
end


function decode(ring::Ring, mx::Array{BigInt, 1}, slots::Int, logp::Int, logq::Int)

    vals = Array{Complex{Float64}}(undef, slots)

    q = ring.qpows[logq+1]
    gap = Nh ÷ slots

    for i in 0:slots-1
        idx = i * gap
        tmp = mod(mx[idx+1], q)

        # Brings `tmp` to the range [-q/2, q/2]
        if num_bits(tmp) == logq
            #println((tmp, logq, q, num_bits(tmp)))
            tmp -= q
        end

        vr = integer_to_float(Float64, tmp, logp)

        tmp = mod(mx[idx + Nh + 1], q)

        # Brings `tmp` to the range [-q/2, q/2]
        if num_bits(tmp) == logq
            #println((tmp, logq, q, num_bits(tmp)))
            tmp -= q
        end

        vi = integer_to_float(Float64, tmp, logp)

        vals[i+1] = vr + im * vi
    end
    EMB(ring, vals, slots)
end


function add(ring::Ring, p1::Array{BigInt, 1}, p2::Array{BigInt, 1}, modulus::BigInt)
    AddMod.(p1, p2, modulus)
end


function sub(ring::Ring, p1::Array{BigInt, 1}, p2::Array{BigInt, 1}, modulus::BigInt)
    AddMod.(p1, -p2, modulus)
end


function multDNTT(ring::Ring, ra::Array{UInt64, 1}, rb::Array{UInt64, 1}, np::Int, q::BigInt)
    multDNTT(ring.multiplier, ra, rb, np, q)
end


function addNTT(ring::Ring, ra::Array{UInt64, 1}, rb::Array{UInt64, 1}, np::Int)
    addNTT(ring.multiplier, ra, rb, np);
end


function multByMonomial(ring::Ring, p::Array{BigInt, 1}, monomialDeg::Int)
    shift = monomialDeg % M
    pp = Polynomial(p, true)
    pp_shifted = shift_polynomial(pp, shift)
    pp_shifted.coeffs
end


function leftRotate(ring::Ring, p::Array{BigInt, 1}, r::Int)
    # TODO: what's the difference from multByMonomial()?
    res = similar(p)
    pow = ring.rotGroup[r+1]
    for i in 0:N-1
        ipow = i * pow
        shift = ipow % M
        if shift < N
            res[shift+1] = p[i+1]
        else
            res[shift - N + 1] = -p[i+1]
        end
    end
    res
end


function conjugate(ring::Ring, p::Array{BigInt, 1}, q::BigInt, keep_positive::Bool=false)
    res = similar(p)
    res[0+1] = p[0+1]
    for i in 1:N-1
        res[i+1] = keep_positive ? NegMod(p[N - i + 1], q) : -p[N - i + 1]
    end
    res
end


function squareNTT(ring::Ring, ra::Array{UInt64, 1}, np::Int, q::BigInt)
    squareNTT(ring.multiplier, ra, np, q)
end


function Base.mod(ring::Ring, p::Array{BigInt, 1}, modulus::BigInt)
    res = similar(p)
    for i in 0:N-1
        res[i+1] = mod(p[i+1], modulus)
    end
    res
end


function negate(ring::Ring, p::Array{BigInt, 1})
    res = similar(p)
    for i in 0:N-1
        res[i+1] = -p[i+1]
    end
    res
end


function multByConst(ring::Ring, p::Array{BigInt, 1}, cnst::BigInt, modulus::BigInt)
    res = similar(p)
    for i in 0:N-1
        res[i+1] = MulMod(p[i+1], cnst, modulus)
    end
    res
end


function maxBits(f::Array{BigInt, 1}, n::Int)
    maximum(num_bits.(f))
end


function addBootContext!(ring::Ring, logSlots::Int, logp::Int)
    if haskey(ring.bootContextMap, logSlots)
        return
    end

    slots = 1 << logSlots
    dslots = slots << 1
    logk = logSlots >> 1

    k = 1 << logk
    gap = Nh >> logSlots

    rpvec = Array{Array{UInt64, 1}}(undef, slots)
    rpvecInv = Array{Array{UInt64, 1}}(undef, slots)

    bndvec = Array{Int}(undef, slots)
    bndvecInv = Array{Int}(undef, slots)

    pvec = Array{BigInt}(undef, N)
    pvec .= 0
    pvals = Array{Complex{Float64}}(undef, dslots)

    c = 0.25 / pi

    if logSlots < logNh
        dgap = gap >> 1
        for ki in 0:k:slots-1
            for pos in ki:ki+k-1
                for i in 0:slots-pos-1
                    deg = ((M - ring.rotGroup[i + pos + 1]) * i * gap) % M
                    pvals[i+1] = ring.ksiPows[deg + 1]
                    pvals[i + slots + 1] = pvals[i+1] * im
                end
                for i in slots-pos:slots-1
                    deg = ((M - ring.rotGroup[i + pos - slots + 1]) * i * gap) % M
                    pvals[i + 1] = ring.ksiPows[deg + 1]
                    pvals[i + slots + 1] = pvals[i + 1] * im
                end
                # TODO: check that this is equivalent to rightRotateAndEqual(pvals, dslots, ki)
                pvals = circshift(pvals, ki)
                pvals = EMBInv(ring, pvals, dslots)
                for i in 0:dslots-1
                    jdx = Nh + i * dgap
                    idx = i * dgap
                    pvec[idx + 1] = float_to_integer(real(pvals[i + 1]), logp)
                    pvec[jdx + 1] = float_to_integer(imag(pvals[i + 1]), logp)
                end
                bndvec[pos + 1] = maxBits(pvec, N)
                np = cld(bndvec[pos + 1] + logQ + 2 * logN + 2, 59)
                rpvec[pos + 1] = CRT(ring, pvec, np)
                for i in 0:N-1
                    pvec[i + 1] = zero(BigInt)
                end
            end
        end

        for i in 0:slots-1
            pvals[i + 1] = 0.0
            pvals[i + slots + 1] = -c * im
        end
        pvals = EMBInv(ring, pvals, dslots)
        for i in 0:dslots-1
            idx = i * dgap
            jdx = Nh + i * dgap
            pvec[idx + 1] = float_to_integer(real(pvals[i+1]), logp)
            pvec[jdx + 1] = float_to_integer(imag(pvals[i+1]), logp)
        end
        bnd1 = maxBits(pvec, N)
        np = cld(bnd1 + logQ + 2 * logN + 2, 59)
        rp1 = CRT(ring, pvec, np)
        for i in 0:N-1
            pvec[i+1] = zero(BigInt)
        end

        for i in 0:slots-1
            pvals[i + 1] = c
            pvals[i + slots + 1] = 0
        end

        pvals = EMBInv(ring, pvals, dslots)
        for i in 0:dslots-1
            idx = i * dgap
            jdx = Nh + i * dgap
            pvec[idx+1] = float_to_integer(real(pvals[i+1]), logp)
            pvec[jdx+1] = float_to_integer(imag(pvals[i+1]), logp)
        end
        bnd2 = maxBits(pvec, N)
        np = cld(bnd2 + logQ + 2 * logN + 2, 59)
        rp2 = CRT(ring, pvec, np)
        for i in 0:N-1
            pvec[i+1] = zero(BigInt)
        end
    else
        # TODO: need to test this branch
        for ki in 0:k:slots-1
            for pos in ki:ki+k-1
                for i in 0:slots-pos-1
                    deg = ((M - ring.rotGroup[i + pos + 1]) * i * gap) % M
                    pvals[i+1] = ring.ksiPows[deg+1]
                end
                for i in slots-pos:slots-1
                    deg = ((M - ring.rotGroup[i + pos - slots + 1]) * i * gap) % M
                    pvals[i+1] = ring.ksiPows[deg]
                end
                # TODO: check that this is equivalent to rightRotateAndEqual(pvals, slots, ki)
                # TODO: in the original it was `slots`, but length of `pvals` is `dslots` - bug?
                pvals = vcat(circshift(pvals[1:slots], ki), pvals[slots+1:end])
                pvals = EMBInv(ring, pvals, slots)
                for i in 0:slots-1
                    idx = i * gap
                    jdx = Nh + i * gap
                    pvec[idx+1] = float_to_integer(real(pvals[i+1]), logp)
                    pvec[jdx+1] = float_to_integer(imag(pvals[i+1]), logp)
                end
                bndvec[pos+1] = maxBits(pvec, N)
                np = cld(bndvec[pos+1] + logQ + 2 * logN + 2, 59)
                rpvec[pos+1] = CRT(ring, pvec, np)
                for i in 0:N-1
                    pvec[i+1] = zero(BigInt)
                end
            end
        end

    end

    for ki in 0:k:slots-1
        for pos in ki:ki+k-1
            for i in 0:slots-pos-1
                deg = (ring.rotGroup[i+1] * (i + pos) * gap) % M
                pvals[i+1] = ring.ksiPows[deg+1]
            end
            for i in slots-pos:slots-1
                deg = (ring.rotGroup[i+1] * (i + pos - slots) * gap) % M
                pvals[i+1] = ring.ksiPows[deg+1]
            end
            # TODO: check that this is equivalent to rightRotateAndEqual(pvals, slots, ki)
            # TODO: in the original it was `slots`, but length of `pvals` is `dslots` - bug?
            pvals = vcat(circshift(pvals[1:slots], ki), pvals[slots+1:end])
            pvals = EMBInv(ring, pvals, slots);
            for i in 0:slots-1
                idx = i * gap
                jdx = Nh + i * gap
                pvec[idx+1] = float_to_integer(real(pvals[i+1]), logp)
                pvec[jdx+1] = float_to_integer(imag(pvals[i+1]), logp)
            end
            bndvecInv[pos+1] = maxBits(pvec, N)
            np = cld(bndvecInv[pos+1] + logQ + 2 * logN + 2, 59)
            rpvecInv[pos+1] = CRT(ring, pvec, np)
            for i in 0:N-1
                pvec[i+1] = zero(BigInt)
            end
        end
    end

    ring.bootContextMap[logSlots] = BootContext(
        rpvec, rpvecInv, rp1, rp2, bndvec, bndvecInv, bnd1, bnd2, logp)
end
