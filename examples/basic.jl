using HEAAN


function coinciding_bits(x::Float64, y::Float64)
    log_diff = floor(Int, -log2(abs(x - y)))
    log_diff >= 0 ? log_diff : 0
end


coinciding_bits(x::Array{Complex{Float64}}, y::Array{Complex{Float64}}) =
    vcat(coinciding_bits.(real.(x), real.(y)), coinciding_bits.(imag.(x), imag.(y)))


function make_compatible(cipher1, cipher2)
    p1 = cipher1.log_precision
    p2 = cipher2.log_precision
    if p1 < p2
        cipher2 = rescale_by(cipher2, p2 - p1)
    elseif p2 < p1
        cipher1 = rescale_by(cipher1, p1 - p2)
    end

    q1 = cipher1.log_cap
    q2 = cipher2.log_cap

    if q1 < q2
        cipher2 = mod_down_by(cipher2, q2 - q1)
    elseif q2 < q1
        cipher1 = mod_down_by(cipher1, q1 - q2)
    end

    cipher1, cipher2
end


function test()
    # Vector size
    # The maximum vector size is polynomial_length / 2
    n = 2^6

    # Initial precision for ciphertexts (that is, the absolute precision is 1/2^30)
    log_precision = 30

    # Initial precision cap for ciphertexts
    # Gets consumed during e.g. multiplication
    log_cap = 200

    # Polynomial length and full modulus determine the security of the scheme
    params = Params(log_polynomial_length=8, log_lo_modulus=300)

    rng = HEAAN.MyRNG(12345)

    secret_key = SecretKey(rng, params)

    enc_key = EncryptionKey(rng, secret_key)
    mul_key = MultiplicationKey(rng, secret_key)
    conj_key = ConjugationKey(rng, secret_key)

    # Rotation keys must be created for each required shift value
    lr_key = LeftRotationKey(rng, secret_key, 16)

    v1 = HEAAN.randomComplexArray(rng, n) # randn(rng, n) + im * randn(rng, n)
    v2 = HEAAN.randomComplexArray(rng, n) # randn(rng, n) + im * randn(rng, n)

    # Reference calculation
    sigmoid_ref(v) = exp(v) / (1 + exp(v))
    ref = circshift(sigmoid_ref.(v1) .* conj.(v2), -16)

    # Encrypt the initial vectors
    c1 = encrypt(rng, enc_key, v1, log_precision, log_cap)
    c2 = encrypt(rng, enc_key, v2, log_precision, log_cap)

    println("Before: precision=$(c1.log_precision), cap=$(c2.log_cap)")

    # Intermediate results
    t1 = sigmoid(mul_key, c1, log_precision, 7)
    t2 = conj(conj_key, c2)
    println("t1: $((t1.log_precision, t1.log_cap)), t2: $((t2.log_precision, t2.log_cap))")

    # The intermediate results have differing `log_cap`
    # (and may also have different `log_precision`, depending on applied functions).
    # In order to perform mul() or add() they must be made equal
    # (technically, `mul()` can handle different `log_precision`, but not `log_cap`)
    t1, t2 = make_compatible(t1, t2)

    println("t1: $((t1.log_precision, t1.log_cap)), t2: $((t2.log_precision, t2.log_cap))")

    # Finish the computation and decrypt
    cres = circshift(lr_key, mul(mul_key, t1, t2), -16)
    res = decrypt(secret_key, cres)

    println("After: precision=$(cres.log_precision), cap=$(cres.log_cap)")

    bits = coinciding_bits(ref, res)
    println("Coinciding bits: $(minimum(bits)) to $(maximum(bits))")
end


test()
