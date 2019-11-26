const ENCRYPTION = 0
const MULTIPLICATION = 1
const CONJUGATION = 2


struct Scheme

    ring :: Ring

    isSerialized :: Bool

    # contain Encryption, Multiplication and Conjugation keys, if generated
    keyMap :: Dict{Int, Key}

    # contain left rotation keys, if generated
    leftRotKeyMap :: Dict{Int, Key}

    # contain Encryption, Multiplication and Conjugation keys, if generated
    serKeyMap :: Dict{Int, String}

    # contain left rotation keys, if generated
    serLeftRotKeyMap :: Dict{Int, String}

    function Scheme(rng::MyRNG, secretKey::SecretKey, ring::Ring, isSerialized::Bool = false)
        enc_key = addEncKey(rng, secretKey, ring)
        mult_key = addMultKey(rng, secretKey, ring)

        keyMap = Dict(ENCRYPTION => enc_key, MULTIPLICATION => mult_key)

        new(ring, isSerialized, keyMap, Dict{Int, Key}())
    end

end


function addEncKey(rng::MyRNG, secretKey::SecretKey, ring::Ring)

    np = cld(1 + logQQ + logN + 2, 59)
    ax = sampleUniform2(rng, logQQ)
    bx = mult(ring, secretKey.sx, ax, np, QQ)
    bx = subFromGaussAndEqual(rng, ring, bx, QQ)

    key = Key()
    key.rax .= CRT(ring, ax, nprimes, logQQ)
    key.rbx .= CRT(ring, bx, nprimes, logQQ)

    key
end


function addMultKey(rng::MyRNG, secretKey::SecretKey, ring::Ring)

    np = cld(1 + logQQ + logN + 2, 59)
    ax = sampleUniform2(rng, logQQ)

    bx = mult(ring, secretKey.sx, ax, np, QQ)
    bx = subFromGaussAndEqual(rng, ring, bx, QQ)

    np = cld(2 + logN + 2, 59)
    sxsx = mult(ring, secretKey.sx, secretKey.sx, np, Q)
    leftShiftAndEqual!(sxsx, logQ, QQ)
    addAndEqual!(bx, sxsx, QQ)

    key = Key()
    key.rax .= CRT(ring, ax, nprimes, logQQ)
    key.rbx .= CRT(ring, bx, nprimes, logQQ)

    key
end


function encryptSingle(rng::MyRNG, scheme::Scheme, val::Complex{Float64}, logp::Int, logq::Int)
    plain = encodeSingle(val, logp, logq)
    encryptMsg(rng, scheme, plain)
end


function encodeSingle(val::Float64, logp::Int, logq::Int)
    p = Plaintext(logp, logq, 1)
    p.mx[0+1] = float_to_integer(val, logp + logQ, logq + logQ)
    p
end


function encodeSingle(val::Complex{Float64}, logp::Int, logq::Int)
    p = Plaintext(logp, logq, 1)
    p.mx .= 0 # TODO: it seems that in the C++ version this happens automatically
    p.mx[0+1] = float_to_integer(real(val), logp + logQ, logq + logQ)
    p.mx[Nh+1] = float_to_integer(imag(val), logp + logQ, logq + logQ)
    p
end


function encrypt(
        rng::MyRNG, scheme::Scheme, vals::Array{Complex{Float64}, 1}, n::Int, logp::Int, logq::Int)
    plain = encode(scheme, vals, n, logp, logq)
    encryptMsg(rng, scheme, plain)
end


function encode(scheme::Scheme, vals::Array{Complex{Float64}, 1}, n::Int, logp::Int, logq::Int)
    p = Plaintext(logp, logq, n)
    p.mx .= 0
    # TODO: note that it only fills certain indices, not the whole array
    encode(scheme.ring, p.mx, vals, n, logp + logQ, logq + logQ)
    p
end


function encryptMsg(rng::MyRNG, scheme::Scheme, plain::Plaintext)

    qQ = scheme.ring.qpows[plain.logq + logQ + 1]

    vx = sampleZO(rng, scheme.ring)

    key = scheme.keyMap[ENCRYPTION]

    np = cld(1 + logQQ + logN + 2, 59)

    ax = multNTT(scheme.ring.multiplier, vx, key.rax, np, qQ)
    ax = addGaussAndEqual(rng, scheme.ring, ax, qQ)

    bx = multNTT(scheme.ring.multiplier, vx, key.rbx, np, qQ)
    bx = addGaussAndEqual(rng, scheme.ring, bx, qQ)

    addAndEqual!(bx, plain.mx, qQ)
    rightShiftAndEqual!(ax, logQ)
    rightShiftAndEqual!(bx, logQ)

    cipher = Ciphertext(plain.logp, plain.logq, plain.n)
    cipher.ax .= ax
    cipher.bx .= bx

    cipher
end


function decryptSingle(scheme::Scheme, secretKey::SecretKey, cipher::Ciphertext)
    plain = decryptMsg(scheme, secretKey, cipher)
    decodeSingle(scheme.ring, plain)
end

function decrypt(scheme::Scheme, secretKey::SecretKey, cipher::Ciphertext)
    plain = decryptMsg(scheme, secretKey, cipher)
    decode(scheme, plain)
end


function decryptMsg(scheme::Scheme, secretKey::SecretKey, cipher::Ciphertext)
    q = scheme.ring.qpows[cipher.logq+1]
    plain = Plaintext(cipher.logp, cipher.logq, cipher.n)
    np = cld(1 + cipher.logq + logN + 2, 59)
    plain.mx .= mult(scheme.ring, cipher.ax, secretKey.sx, np, q)
    addAndEqual!(plain.mx, cipher.bx, q)
    plain
end


function decodeSingle(ring::Ring, plain::Plaintext)
    q = ring.qpows[plain.logq+1]

    tmp = plain.mx[0+1] % q
    if num_bits(tmp) == plain.logq
        tmp -= q
    end
    r = integer_to_float(Float64, tmp, plain.logp)

    tmp = plain.mx[Nh+1] % q
    if num_bits(tmp) == plain.logq
        tmp -= q
    end
    i = integer_to_float(Float64, tmp, plain.logp)

    r + im * i
end


function decode(scheme::Scheme, plain::Plaintext)
    decode(scheme.ring, plain.mx, plain.n, plain.logp, plain.logq)
end


function add(scheme::Scheme, cipher1::Ciphertext, cipher2::Ciphertext)
    q = scheme.ring.qpows[cipher1.logq+1]
    cipher_res = Ciphertext(cipher1.logp, cipher1.logq, cipher1.n)
    cipher_res.ax .= add(scheme.ring, cipher1.ax, cipher2.ax, q)
    cipher_res.bx .= add(scheme.ring, cipher1.bx, cipher2.bx, q)
    cipher_res
end


function sub(scheme::Scheme, cipher1::Ciphertext, cipher2::Ciphertext)
    q = scheme.ring.qpows[cipher1.logq+1]
    cipher_res = Ciphertext(cipher1.logp, cipher1.logq, cipher1.n)
    cipher_res.ax .= sub(scheme.ring, cipher1.ax, cipher2.ax, q)
    cipher_res.bx .= sub(scheme.ring, cipher1.bx, cipher2.bx, q)
    cipher_res
end


function mult(scheme::Scheme, cipher1::Ciphertext, cipher2::Ciphertext)
    res = Ciphertext(cipher1.logp + cipher2.logp, cipher1.logq, cipher1.n)

    @assert cipher1.logq == cipher2.logq

    ring = scheme.ring

    q = ring.qpows[cipher1.logq+1]
    qQ = ring.qpows[cipher1.logq + logQ + 1]

    np = cld(2 + cipher1.logq + cipher2.logq + logN + 2, 59)

    ra1 = CRT(ring, cipher1.ax, np, cipher1.logq)
    rb1 = CRT(ring, cipher1.bx, np, cipher1.logq)
    ra2 = CRT(ring, cipher2.ax, np, cipher2.logq)
    rb2 = CRT(ring, cipher2.bx, np, cipher2.logq)

    # TODO: shouldn't it be logq1 + logq2? Or need to check that logq1==logq2
    axax = multDNTT(ring, ra1, ra2, np, q)
    bxbx = multDNTT(ring, rb1, rb2, np, q)

    ra1 = addNTT(ring, ra1, rb1, np)
    ra2 = addNTT(ring, ra2, rb2, np)

    axbx = multDNTT(ring, ra1, ra2, np, q)

    key = scheme.keyMap[MULTIPLICATION]

    np = cld(cipher1.logq + logQQ + logN + 2, 59)
    raa = CRT(ring, axax, np, cipher1.logq) # TODO: see the note above
    res.ax .= multDNTT(ring, raa, key.rax, np, qQ)
    res.bx .= multDNTT(ring, raa, key.rbx, np, qQ)

    rightShiftAndEqual!(res.ax, logQ)
    rightShiftAndEqual!(res.bx, logQ)

    res.ax .= add(ring, res.ax, axbx, q)
    res.ax .= sub(ring, res.ax, bxbx, q)
    res.ax .= sub(ring, res.ax, axax, q)
    res.bx .= add(ring, res.bx, bxbx, q)

    res
end


function imult(scheme::Scheme, cipher::Ciphertext)
    res = Ciphertext(cipher.logp, cipher.logq, cipher.n)
    res.ax .= multByMonomial(scheme.ring, cipher.ax, Nh)
    res.bx .= multByMonomial(scheme.ring, cipher.bx, Nh)
    res
end


function multByPolyNTT(scheme::Scheme, cipher::Ciphertext, rpoly::Array{UInt64, 1}, bnd::Int, logp::Int)
    ring = scheme.ring
    q = ring.qpows[cipher.logq+1]

    res = Ciphertext(cipher.logp + logp, cipher.logq, cipher.n)

    np = cld(cipher.logq + bnd + logN + 2, 59)
    res.ax .= multNTT(ring.multiplier, cipher.ax, rpoly, np, q)
    res.bx .= multNTT(ring.multiplier, cipher.bx, rpoly, np, q)

    res
end


function addLeftRotKey!(rng::MyRNG, scheme::Scheme, secretKey::SecretKey, r::Int)

    ring = scheme.ring

    np = cld(1 + logQQ + logN + 2, 59)
    ax = sampleUniform2(rng, logQQ)
    bx = mult(ring, secretKey.sx, ax, np, QQ)
    bx = subFromGaussAndEqual(rng, ring, bx, QQ)

    spow = leftRotate(ring, secretKey.sx, r, QQ)
    leftShiftAndEqual!(spow, logQ, QQ)
    bx = add(ring, bx, spow, QQ)

    key = Key()
    key.rax .= CRT(ring, ax, nprimes, logQQ)
    key.rbx .= CRT(ring, bx, nprimes, logQQ)

    push!(scheme.leftRotKeyMap, r => key)
end


function addLeftRotKeys!(rng::MyRNG, scheme::Scheme, secretKey::SecretKey)
    for i in 0:logN-2
        idx = 1 << i
        if !haskey(scheme.leftRotKeyMap, idx)
            addLeftRotKey!(rng, scheme, secretKey, idx)
        end
    end
end


function leftRotateFast(scheme::Scheme, cipher::Ciphertext, r::Int)

    ring = scheme.ring

    q = ring.qpows[cipher.logq+1]
    qQ = ring.qpows[cipher.logq + logQ + 1]

    cipher_res = Ciphertext(cipher.logp, cipher.logq, cipher.n)

    check_range(cipher.ax, cipher.logq)
    check_range(cipher.bx, cipher.logq)

    bxrot = leftRotate(ring, cipher.bx, r, one(BigInt) << cipher.logq)
    axrot = leftRotate(ring, cipher.ax, r, one(BigInt) << cipher.logq)

    check_range(bxrot, cipher.logq)
    check_range(axrot, cipher.logq)

    key = scheme.leftRotKeyMap[r]
    np = cld(cipher.logq + logQQ + logN + 2, 59)

    rarot = CRT(ring, axrot, np, cipher.logq)

    cipher_res.ax .= multDNTT(ring, rarot, key.rax, np, qQ)
    cipher_res.bx .= multDNTT(ring, rarot, key.rbx, np, qQ)

    rightShiftAndEqual!(cipher_res.ax, logQ)
    rightShiftAndEqual!(cipher_res.bx, logQ)

    # TODO: Some of `bxrot` elements are >q here (in the original as well)
    # Need to check and see where it is coming from.
    cipher_res.bx .= add(ring, cipher_res.bx, bxrot, q)

    check_range(cipher_res.ax, cipher_res.logq)
    check_range(cipher_res.bx, cipher_res.logq)

    cipher_res
end


function addConjKey!(rng::MyRNG, scheme::Scheme, secretKey::SecretKey)
    ring = scheme.ring

    np = cld(1 + logQQ + logN + 2, 59)
    ax = sampleUniform2(rng, logQQ)
    check_range(ax, logQQ)
    bx = mult(ring, secretKey.sx, ax, np, QQ)
    check_range(bx, logQQ)
    bx = subFromGaussAndEqual(rng, ring, bx, QQ)
    check_range(bx, logQQ)

    sxconj = conjugate(ring, secretKey.sx, QQ, true)
    check_range(sxconj, logQQ)
    leftShiftAndEqual!(sxconj, logQ, QQ)
    check_range(sxconj, logQQ)
    bx = add(ring, bx, sxconj, QQ)
    check_range(bx, logQQ)

    key = Key()
    key.rax .= CRT(ring, ax, nprimes, logQQ)
    key.rbx .= CRT(ring, bx, nprimes, logQQ)

    scheme.keyMap[CONJUGATION] = key
end


function conjugate(scheme::Scheme, cipher::Ciphertext)

    ring = scheme.ring

    q = ring.qpows[cipher.logq+1]
    qQ = ring.qpows[cipher.logq + logQ + 1]

    cipher_res = Ciphertext(cipher.logp, cipher.logq, cipher.n)

    check_range(cipher.ax, cipher.logq)
    check_range(cipher.bx, cipher.logq)

    bxconj = conjugate(ring, cipher.bx, one(BigInt) << cipher.logq, true)

    # TODO: setting `true` here loses compatibility with C++ HEAAN,
    # but keeps the precision the same.
    axconj = conjugate(ring, cipher.ax, one(BigInt) << cipher.logq, true)

    key = scheme.keyMap[CONJUGATION]

    np = cld(cipher.logq + logQQ + logN + 2, 59)
    raconj = CRT(ring, axconj, np, cipher.logq)
    cipher_res.ax .= multDNTT(ring, raconj, key.rax, np, qQ)
    cipher_res.bx .= multDNTT(ring, raconj, key.rbx, np, qQ)

    check_range(cipher_res.ax, cipher.logq + logQ)
    check_range(cipher_res.bx, cipher.logq + logQ)

    rightShiftAndEqual!(cipher_res.ax, logQ)
    rightShiftAndEqual!(cipher_res.bx, logQ)

    check_range(cipher_res.ax, cipher.logq)
    check_range(cipher_res.bx, cipher.logq)

    cipher_res.bx .= add(ring, cipher_res.bx, bxconj, q)
    check_range(cipher_res.bx, cipher.logq)

    cipher_res
end


function square(scheme::Scheme, cipher::Ciphertext)

    # TODO: can be joined with mult()?

    ring = scheme.ring

    res = Ciphertext(cipher.logp * 2, cipher.logq, cipher.n)

    q = ring.qpows[cipher.logq + 1]
    qQ = ring.qpows[cipher.logq + logQ + 1]

    np = cld(2 * cipher.logq + logN + 2, 59)

    ra = CRT(ring, cipher.ax, np, cipher.logq)
    rb = CRT(ring, cipher.bx, np, cipher.logq)

    bxbx = squareNTT(ring, rb, np, q)
    axax = squareNTT(ring, ra, np, q)
    axbx = multDNTT(ring, ra, rb, np, q)
    axbx = add(ring, axbx, axbx, q)

    key = scheme.keyMap[MULTIPLICATION]

    np = cld(cipher.logq + logQQ + logN + 2, 59)
    raa = CRT(ring, axax, np, cipher.logq)
    res.ax .= multDNTT(ring, raa, key.rax, np, qQ)
    res.bx .= multDNTT(ring, raa, key.rbx, np, qQ)

    rightShiftAndEqual!(res.ax, logQ)
    rightShiftAndEqual!(res.bx, logQ)

    res.ax .= add(ring, res.ax, axbx, q)
    res.bx .= add(ring, res.bx, bxbx, q)

    res
end


function reScaleBy(scheme::Scheme, cipher::Ciphertext, dlogq::Int)

    res = Ciphertext(cipher.logp - dlogq, cipher.logq - dlogq, cipher.n)

    res.ax .= cipher.ax
    res.bx .= cipher.bx
    rightShiftAndEqual!(res.ax, dlogq)
    rightShiftAndEqual!(res.bx, dlogq)

    res
end


function modDownBy(scheme::Scheme, cipher::Ciphertext, dlogq::Int)
    ring = scheme.ring
    q = ring.qpows[cipher.logq - dlogq + 1]

    res = Ciphertext(cipher.logp, cipher.logq - dlogq, cipher.n)
    res.ax .= mod(ring, cipher.ax, q)
    res.bx .= mod(ring, cipher.bx, q)

    res
end


function modDownTo(scheme::Scheme, cipher::Ciphertext, logq::Int)
    ring = scheme.ring
    q = ring.qpows[logq + 1]

    res = Ciphertext(cipher.logp, logq, cipher.n)
    res.ax .= mod(ring, cipher.ax, q)
    res.bx .= mod(ring, cipher.bx, q)

    res
end


function negate(scheme::Scheme, cipher::Ciphertext)
    ring = scheme.ring
    res = Ciphertext(cipher.logp, cipher.logq, cipher.n)
    res.ax .= negate(ring, cipher.ax)
    res.bx .= negate(ring, cipher.bx)
    res
end


function addConst(scheme::Scheme, cipher::Ciphertext, cnst::Union{BigFloat, Float64}, logp::Int)
    ring = scheme.ring
    q = ring.qpows[cipher.logq + 1]
    res = copy(cipher)
    cnstZZ = logp < 0 ? float_to_integer(cnst, cipher.logp, cipher.logq) : float_to_integer(cnst, logp, cipher.logq)
    res.bx[0+1] = AddMod(cipher.bx[0+1], cnstZZ, q)
    res
end


function multByConst(scheme::Scheme, cipher::Ciphertext, cnst::Union{BigFloat, Float64}, logp::Int)
    ring = scheme.ring
    q = ring.qpows[cipher.logq + 1]
    cnstZZ = float_to_integer(cnst, logp, cipher.logq)

    res = Ciphertext(cipher.logp + logp, cipher.logq, cipher.n)
    res.ax .= multByConst(ring, cipher.ax, cnstZZ, q)
    res.bx .= multByConst(ring, cipher.bx, cnstZZ, q)

    res
end


function multByConst(scheme::Scheme, cipher::Ciphertext, cnst::Complex{Float64}, logp::Int)
    ring = scheme.ring
    q = ring.qpows[cipher.logq + 1]
    cnstZZReal = float_to_integer(real(cnst), logp, cipher.logq)
    cnstZZImag = float_to_integer(imag(cnst), logp, cipher.logq)

    tmp = Ciphertext(cipher.logp, cipher.logq, cipher.n)
    tmp.ax .= multByConst(ring, cipher.ax, cnstZZImag, q)
    tmp.bx .= multByConst(ring, cipher.bx, cnstZZImag, q)
    tmp.ax .= multByMonomial(ring, tmp.ax, N ÷ 2)
    tmp.bx .= multByMonomial(ring, tmp.bx, N ÷ 2)

    res = Ciphertext(cipher.logp + logp, cipher.logq, cipher.n)
    res.ax .= multByConst(ring, cipher.ax, cnstZZReal, q)
    res.bx .= multByConst(ring, cipher.bx, cnstZZReal, q)
    addAndEqual!(res.ax, tmp.ax, QQ)
    addAndEqual!(res.bx, tmp.bx, QQ)

    res
end


function multByConstVec(scheme::Scheme, cipher::Ciphertext, cnstVec::Array{Complex{Float64}, 1}, logp::Int)
    slots = cipher.n
    cnstPoly = Array{BigInt}(undef, N)
    encode(ring, cnstPoly, cnstVec, slots, logp)
    multByPoly(scheme, cipher, cnstPoly, logp)
end


function multByPoly(scheme::Scheme, cipher::Ciphertext, poly::Array{BigInt, 1}, logp::Int)
    ring = scheme.ring
    q = ring.qpows[cipher.logq + 1]
    res = Ciphertext(cipher.logp + logp, cipher.logq, cipher.n)

    bnd = maxBits(poly, N)
    np = cld(cipher.logq + bnd + logN + 2, pbnd)
    rpoly = CRT(ring, poly, np, cipher.logq)
    res.ax .= multNTT(ring.multiplier, cipher.ax, rpoly, np, q)
    res.bx .= multNTT(ring.multiplier, cipher.bx, rpoly, np, q)

    res
end


function addBootKey!(rng::MyRNG, scheme::Scheme, secretKey::SecretKey, logl::Int, logp::Int)

    ring = scheme.ring

    addBootContext!(ring, logl, logp)
    addConjKey!(rng, scheme, secretKey)
    addLeftRotKeys!(rng, scheme, secretKey)

    loglh = logl ÷ 2
    k = 1 << loglh
    m = 1 << (logl - loglh)

    # TODO: isn't it what addLeftRotKeys!() does?
    for i in 1:k-1
        if !haskey(scheme.leftRotKeyMap, i)
            addLeftRotKey!(rng, scheme, secretKey, i)
        end
    end

    for i in 1:m-1
        idx = i * k
        if !haskey(scheme.leftRotKeyMap, idx)
            addLeftRotKey!(rng, scheme, secretKey, idx)
        end
    end
end


function normalize(scheme::Scheme, cipher::Ciphertext)
    # Since we want all numbers to be positive, this does nothing
    return cipher

    q = scheme.ring.qpows[cipher.logq+1]
    new_cipher = Ciphertext(cipher.logp, cipher.logq, cipher.n)
    new_cipher.ax .= cipher.ax
    new_cipher.bx .= cipher.bx
    for i in 0:N-1
        if num_bits(new_cipher.ax[i+1]) == new_cipher.logq
            new_cipher.ax[i+1] -= q
        end
        if num_bits(new_cipher.bx[i+1]) == new_cipher.logq
            new_cipher.bx[i+1] -= q
        end
    end
    new_cipher
end


function divByPo2(scheme::Scheme, cipher::Ciphertext, bits::Int)
    new_cipher = Ciphertext(cipher.logp, cipher.logq - bits, cipher.n)
    new_cipher.ax .= cipher.ax
    new_cipher.bx .= cipher.bx
    rightShiftAndEqual!(new_cipher.ax, bits)
    rightShiftAndEqual!(new_cipher.bx, bits)
    new_cipher
end


function coeffToSlot(scheme::Scheme, cipher::Ciphertext)

    ring = scheme.ring

    slots = cipher.n
    logSlots = floor(Int, log2(slots)) # TODO: check it's actually `floor` and not `round`
    logk = logSlots ÷ 2
    k = 1 << logk

    rotvec = Array{Ciphertext}(undef, k)
    rotvec[0+1] = copy(cipher)

    for j in 0:k-2
        rotvec[j+1+1] = leftRotateFast(scheme, rotvec[0+1], j + 1)
    end

    bootContext = ring.bootContextMap[logSlots]

    tmpvec = Array{Ciphertext}(undef, k)

    for j in 0:k-1
        tmpvec[j+1] = multByPolyNTT(
            scheme, rotvec[j+1], bootContext.rpvec[j+1], bootContext.bndvec[j+1], bootContext.logp)
    end

    for j in 1:k-1
        tmpvec[0+1] = add(scheme, tmpvec[0+1], tmpvec[j+1])
    end

    cipher = copy(tmpvec[0+1])

    for ki in k:k:slots-1
        for j in 0:k-1
            tmpvec[j+1] = multByPolyNTT(
                scheme, rotvec[j+1], bootContext.rpvec[j+ki+1], bootContext.bndvec[j+ki+1],
                bootContext.logp)
        end
        for j in 1:k-1
            tmpvec[0+1] = add(scheme, tmpvec[0+1], tmpvec[j+1])
        end
        tmpvec[0+1] = leftRotateFast(scheme, tmpvec[0+1], ki)
        cipher = add(scheme, cipher, tmpvec[0+1])
    end
    reScaleBy(scheme, cipher, bootContext.logp)
end


function slotToCoeff(scheme::Scheme, cipher::Ciphertext)
    ring = scheme.ring

    slots = cipher.n
    logSlots = floor(Int, log2(slots)) # TODO: check hat it is actually floor and not round
    logk = logSlots ÷ 2
    k = 1 << logk

    rotvec = Array{Ciphertext}(undef, k)
    rotvec[0+1] = copy(cipher)

    for j in 0:k-1-1
        rotvec[j+1+1] = leftRotateFast(scheme, rotvec[0+1], j + 1)
    end

    bootContext = ring.bootContextMap[logSlots]

    tmpvec = Array{Ciphertext}(undef, k)

    for j in 0:k-1
        tmpvec[j+1] = multByPolyNTT(
            scheme, rotvec[j+1], bootContext.rpvecInv[j+1],
            bootContext.bndvecInv[j+1], bootContext.logp)
    end

    for j in 1:k-1
        tmpvec[0+1] = add(scheme, tmpvec[0+1], tmpvec[j+1])
    end
    cipher = copy(tmpvec[0+1])

    for ki in k:k:slots-1
        for j in 0:k-1
            tmpvec[j+1] = multByPolyNTT(
                scheme, rotvec[j+1], bootContext.rpvecInv[j+ki+1],
                bootContext.bndvecInv[j+ki+1], bootContext.logp)
        end

        for j in 1:k-1
            tmpvec[0+1] = add(scheme, tmpvec[0+1], tmpvec[j+1])
        end

        tmpvec[0+1] = leftRotateFast(scheme, tmpvec[0+1], ki)
        cipher = add(scheme, cipher, tmpvec[0+1])
    end
    reScaleBy(scheme, cipher, bootContext.logp)
end


function exp2pi(scheme::Scheme, cipher::Ciphertext, logp::Int)
    Pi = BigFloat(pi)

    cipher2 = square(scheme, cipher)
    cipher2 = reScaleBy(scheme, cipher2, logp) # cipher2.logq : logq - logp

    cipher4 = square(scheme, cipher2)
    cipher4 = reScaleBy(scheme, cipher4, logp) # cipher4.logq : logq -2logp
    c = 1/(2 * Pi)
    cipher01 = addConst(scheme, cipher, c, logp) # cipher01.logq : logq

    c = 2*Pi
    cipher01 = multByConst(scheme, cipher01, c, logp)
    cipher01 = reScaleBy(scheme, cipher01, logp) # cipher01.logq : logq - logp

    c = 3/(2*Pi)
    cipher23 = addConst(scheme, cipher, c, logp) # cipher23.logq : logq

    c = 4*Pi*Pi*Pi/3
    cipher23 = multByConst(scheme, cipher23, c, logp)
    cipher23 = reScaleBy(scheme, cipher23, logp) # cipher23.logq : logq - logp

    cipher23 = mult(scheme, cipher23, cipher2)
    cipher23 = reScaleBy(scheme, cipher23, logp) # cipher23.logq : logq - 2logp

    cipher23 = add(scheme, cipher23, cipher01) # cipher23.logq : logq - 2logp

    c = 5/(2*Pi)
    cipher45 = addConst(scheme, cipher, c, logp) # cipher45.logq : logq

    c = 4*Pi*Pi*Pi*Pi*Pi/15
    cipher45 = multByConst(scheme, cipher45, c, logp)
    cipher45 = reScaleBy(scheme, cipher45, logp) # cipher45.logq : logq - logp

    c = 7/(2*Pi)
    cipher = addConst(scheme, cipher, c, logp) # cipher.logq : logq

    c = 8*Pi*Pi*Pi*Pi*Pi*Pi*Pi/315
    cipher = multByConst(scheme, cipher, c, logp)
    cipher = reScaleBy(scheme, cipher, logp) # cipher.logq : logq - logp

    cipher = mult(scheme, cipher, cipher2)
    cipher = reScaleBy(scheme, cipher, logp) # cipher.logq : logq - 2logp

    cipher45 = modDownBy(scheme, cipher45, logp) # cipher45.logq : logq - 2logp
    cipher = add(scheme, cipher, cipher45) # cipher.logq : logq - 2logp

    cipher = mult(scheme, cipher, cipher4)
    cipher = reScaleBy(scheme, cipher, logp) # cipher.logq : logq - 3logp

    cipher23 = modDownBy(scheme, cipher23, logp)
    cipher = add(scheme, cipher, cipher23) # cipher.logq : logq - 3logp

    cipher
end


function evalExp(scheme::Scheme, cipher::Ciphertext, logT::Int, logI::Int=4)
    ring = scheme.ring
    slots = cipher.n
    logSlots = floor(Int, log2(slots)) # TODO: check that it is actually floor and not round
    bootContext = ring.bootContextMap[logSlots]
    if logSlots < logNh
        tmp = conjugate(scheme, cipher)
        check_range(tmp)
        cipher = sub(scheme, cipher, tmp)
        check_range(cipher)
        cipher = divByPo2(scheme, cipher, logT + 1) # bitDown: logT + 1
        check_range(cipher)
        cipher = exp2pi(scheme, cipher, bootContext.logp) # bitDown: logT + 1 + 3(logq + logI)
        check_range(cipher)
        for i in 0:logI+logT-1
            cipher = square(scheme, cipher)
            cipher = reScaleBy(scheme, cipher, bootContext.logp)
        end
        tmp = conjugate(scheme, cipher)
        cipher = sub(scheme, cipher, tmp)
        tmp = multByPolyNTT(scheme, cipher, bootContext.rp1, bootContext.bnd1, bootContext.logp)
        tmprot = leftRotateFast(scheme, tmp, slots)
        tmp = add(scheme, tmp, tmprot)
        cipher = multByPolyNTT(scheme, cipher, bootContext.rp2, bootContext.bnd2, bootContext.logp)
        tmprot = leftRotateFast(scheme, cipher, slots)
        cipher = add(scheme, cipher, tmprot)
        cipher = add(scheme, cipher, tmp)
    else # TODO: check this branch
        tmp = conjugate(scheme, cipher)
        c2 = sub(scheme, cipher, tmp)
        cipher = add(scheme, cipher, tmp)
        cipher = imult(scheme, cipher)
        cihper = divByPo2(scheme, cipher, logT + 1) # cipher bitDown: logT + 1
        c2 = reScaleBy(scheme, c2, logT + 1) # c2 bitDown: logT + 1
        cipher = exp2pi(scheme, cipher, bootContext.logp) # cipher bitDown: logT + 1 + 3(logq + logI)
        c2 = exp2pi(scheme, c2, bootContext.logp) # c2 bitDown: logT + 1 + 3(logq + logI)
        for i in 0:logI+logT-1
            c2 = square(scheme, c2)
            cipher = square(scheme, cipher)
            c2 = reScaleBy(scheme, c2, bootContext.logp)
            cipher = reScaleBy(scheme, cipher, bootContext.logp)
        end
        tmp = conjugate(scheme, c2)
        c2 = sub(scheme, c2, tmp)
        tmp = conjugate(scheme, cipher)
        cipher = sub(scheme, cipher, tmp)
        cipher = imult(scheme, cipher)
        cipher = sub(scheme, c2, cipher)
        c = 0.25/Pi
        cipher = multByConst(scheme, cipher, c, bootContext.logp)
    end
    reScaleBy(scheme, cipher, bootContext.logp + logI)
end
