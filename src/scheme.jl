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

        new(ring, isSerialized, keyMap)
    end

end


function addEncKey(rng::MyRNG, secretKey::SecretKey, ring::Ring)

    np = cld(1 + logQQ + logN + 2, 59)
    ax = sampleUniform2(rng, logQQ)
    bx = mult(ring, secretKey.sx, ax, np, QQ)
    bx = subFromGaussAndEqual(rng, ring, bx, QQ)

    key = Key()
    key.rax .= CRT(ring, ax, nprimes)
    key.rbx .= CRT(ring, bx, nprimes)

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
    key.rax .= CRT(ring, ax, nprimes)
    key.rbx .= CRT(ring, bx, nprimes)

    key
end


function encryptSingle(rng::MyRNG, scheme::Scheme, val::Complex{Float64}, logp::Int, logq::Int)
    plain = encodeSingle(val, logp, logq)
    encryptMsg(rng, scheme, plain)
end


function encodeSingle(val::Float64, logp::Int, logq::Int)
    p = Plaintext(logp, logq, 1)
    p.mx[0+1] = scaleUpToZZ(val, logp + logQ)
    p
end


function encodeSingle(val::Complex{Float64}, logp::Int, logq::Int)
    p = Plaintext(logp, logq, 1)
    p.mx .= 0 # TODO: it seems that in the C++ version this happens automatically
    p.mx[0+1] = scaleUpToZZ(real(val), logp + logQ)
    p.mx[Nh+1] = scaleUpToZZ(imag(val), logp + logQ)
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
    encode(scheme.ring, p.mx, vals, n, logp + logQ)
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
    if NumBits(tmp) == plain.logq
        tmp -= q
    end
    r = scaleDownToReal(tmp, plain.logp)

    tmp = plain.mx[Nh+1] % q
    if NumBits(tmp) == plain.logq
        tmp -= q
    end
    i = scaleDownToReal(tmp, plain.logp)

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
