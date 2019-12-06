using HEAAN
using HEAAN: MyRNG, randomComplexArray, randomCircle, randomComplex, randomReal


function coinciding_bits(x::Float64, y::Float64)
    log_diff = floor(Int, -log2(abs(x - y)))
    log_diff >= 0 ? log_diff : 0
end


coinciding_bits(x::Array{Complex{Float64}}, y::Array{Complex{Float64}}) =
    vcat(coinciding_bits.(real.(x), real.(y)), coinciding_bits.(imag.(x), imag.(y)))


mean(x) = sum(x) / length(x)


std(x) = sqrt(sum((x .- mean(x)).^2) / (length(x) - 1))


function test_approx(x::Array{Complex{Float64}}, y::Array{Complex{Float64}}, expected_min_bits)
    bits = coinciding_bits(x, y)
    min_bits = minimum(bits)
    max_bits = maximum(bits)
    mean_bits = round(mean(bits), digits=2)
    std_bits = round(std(bits), digits=2)
    if min_bits >= expected_min_bits
        @test_result "$min_bits to $max_bits bits ($mean_bits±$std_bits)"
    else
        @test_fail "Minimum coinciding bits: $min_bits, expected $expected_min_bits"
    end
end


@testgroup "Public API" begin


@testcase "encrypt()/decrypt()" begin
    n = 2^6
    log_precision = 30
    log_cap = 100

    #rng = MersenneTwister(12345)
    rng = MyRNG(12345)
    params = Params(log_polynomial_length=8, log_lo_modulus=300)

    secret_key = SecretKey(rng, params)
    enc_key = EncryptionKey(rng, secret_key)

    mvec = randomComplexArray(rng, n) # randn(rng, n) + im * randn(rng, n)
    cipher = encrypt(rng, enc_key, mvec, log_precision, log_cap)
    dvec = decrypt(secret_key, cipher)

    test_approx(mvec, dvec, 23)
end


@testcase "add()" begin
    n = 2^6
    log_precision = 30
    log_cap = 100

    rng = MyRNG(12345)
    params = Params(log_polynomial_length=8, log_lo_modulus=300)

    secret_key = SecretKey(rng, params)

    enc_key = EncryptionKey(rng, secret_key)

    mvec1 = randomComplexArray(rng, n) # randn(rng, n) + im * randn(rng, n)
    mvec2 = randomComplexArray(rng, n) # randn(rng, n) + im * randn(rng, n)

    cipher1 = encrypt(rng, enc_key, mvec1, log_precision, log_cap)
    cipher2 = encrypt(rng, enc_key, mvec2, log_precision, log_cap)

    cipher_res = add(cipher1, cipher2)

    dvec = decrypt(secret_key, cipher_res)

    test_approx(mvec1 .+ mvec2, dvec, 23)
end


@testcase "add_const()" for tp in ["real", "complex"]
    n = 2^6
    log_precision = 30
    log_cap = 100

    rng = MyRNG(12345)
    params = Params(log_polynomial_length=8, log_lo_modulus=300)

    secret_key = SecretKey(rng, params)

    enc_key = EncryptionKey(rng, secret_key)

    mvec = randomComplexArray(rng, n) # randn(rng, n) + im * randn(rng, n)
    mconst = tp == "complex" ? randomComplex(rng) : randomReal(rng) # randn(rng, n) + im * randn(rng, n)

    cipher = encrypt(rng, enc_key, mvec, log_precision, log_cap)

    cipher_res = add_const(cipher, mconst)

    dvec = decrypt(secret_key, cipher_res)

    test_approx(mvec .+ mconst, dvec, 23)
end


@testcase "mul()" begin
    n = 2^6
    log_precision = 30
    log_cap = 100

    rng = MyRNG(12345)
    params = Params(log_polynomial_length=8, log_lo_modulus=300)

    secret_key = SecretKey(rng, params)

    enc_key = EncryptionKey(rng, secret_key)
    mul_key = MultiplicationKey(rng, secret_key)

    mvec1 = randomComplexArray(rng, n) # randn(rng, n) + im * randn(rng, n)
    mvec2 = randomComplexArray(rng, n) # randn(rng, n) + im * randn(rng, n)

    cipher1 = encrypt(rng, enc_key, mvec1, log_precision, log_cap)
    cipher2 = encrypt(rng, enc_key, mvec2, log_precision, log_cap)

    cipher_res = mul(mul_key, cipher1, cipher2)

    dvec = decrypt(secret_key, cipher_res)

    test_approx(mvec1 .* mvec2, dvec, 23)
end


@testcase "mul_by_const()" for tp in ["real", "complex"]
    n = 2^6
    log_precision = 30
    log_cap = 100

    rng = MyRNG(12345)
    params = Params(log_polynomial_length=8, log_lo_modulus=300)

    secret_key = SecretKey(rng, params)

    enc_key = EncryptionKey(rng, secret_key)

    mvec = randomComplexArray(rng, n) # randn(rng, n) + im * randn(rng, n)
    mconst = tp == "complex" ? randomComplex(rng) : randomReal(rng) # randn(rng, n) + im * randn(rng, n)

    cipher = encrypt(rng, enc_key, mvec, log_precision, log_cap)

    cipher_res = mul_by_const(cipher, mconst, log_precision)

    dvec = decrypt(secret_key, cipher_res)

    test_approx(mvec .* mconst, dvec, 23)
end


@testcase "mul_by_const_vec()" begin
    n = 2^6
    log_precision = 30
    log_cap = 300

    rng = MyRNG(12345)
    params = Params(log_polynomial_length=8, log_lo_modulus=300)

    secret_key = SecretKey(rng, params)

    enc_key = EncryptionKey(rng, secret_key)

    mvec = randomComplexArray(rng, n) # randn(rng, n) + im * randn(rng, n)
    cnst = randomComplexArray(rng, n) # randn(rng, n) + im * randn(rng, n)

    cipher = encrypt(rng, enc_key, mvec, log_precision, log_cap)

    cipher_res = mul_by_const_vec(cipher, cnst, log_precision)

    dvec = decrypt(secret_key, cipher_res)

    test_approx(mvec .* cnst, dvec, 23)
end


@testcase "imul()" begin
    n = 2^6
    log_precision = 30
    log_cap = 100

    rng = MyRNG(12345)
    params = Params(log_polynomial_length=8, log_lo_modulus=300)

    secret_key = SecretKey(rng, params)

    enc_key = EncryptionKey(rng, secret_key)

    mvec = randomComplexArray(rng, n) # randn(rng, n) + im * randn(rng, n)

    cipher = encrypt(rng, enc_key, mvec, log_precision, log_cap)

    cipher_res = imul(cipher)

    dvec = decrypt(secret_key, cipher_res)

    test_approx(mvec .* im, dvec, 23)
end


@testcase "circshift()" begin
    n = 2^6
    r = 2^4
    log_precision = 30
    log_cap = 100

    rng = MyRNG(12345)
    params = Params(log_polynomial_length=8, log_lo_modulus=300)

    secret_key = SecretKey(rng, params)

    enc_key = EncryptionKey(rng, secret_key)
    rot_key = LeftRotationKey(rng, secret_key, r)

    mvec = randomComplexArray(rng, n) # randn(rng, n) + im * randn(rng, n)

    cipher = encrypt(rng, enc_key, mvec, log_precision, log_cap)

    cipher_res = circshift(rot_key, cipher, -r)

    dvec = decrypt(secret_key, cipher_res)

    test_approx(circshift(mvec, -r), dvec, 23)
end


@testcase "conj()" begin
    n = 2^6
    log_precision = 30
    log_cap = 300

    rng = MyRNG(12345)
    params = Params(log_polynomial_length=8, log_lo_modulus=300)

    secret_key = SecretKey(rng, params)

    enc_key = EncryptionKey(rng, secret_key)
    conj_key = ConjugationKey(rng, secret_key)

    mvec = randomComplexArray(rng, n) # randn(rng, n) + im * randn(rng, n)

    cipher = encrypt(rng, enc_key, mvec, log_precision, log_cap)

    cipher_res = conj(conj_key, cipher)

    dvec = decrypt(secret_key, cipher_res)

    test_approx(conj.(mvec), dvec, 21)
end


@testcase "power_of_2()" begin
    n = 2^6
    log_degree = 4
    degree = 2^log_degree
    log_precision = 30
    log_cap = 300

    rng = MyRNG(12345)
    params = Params(log_polynomial_length=8, log_lo_modulus=300)

    secret_key = SecretKey(rng, params)

    enc_key = EncryptionKey(rng, secret_key)
    mul_key = MultiplicationKey(rng, secret_key)

    mvec = [randomCircle(rng) for i in 1:n] # randn(rng, n) + im * randn(rng, n)

    # TODO: Temporary replacement for bit-to-bit compatibility with the reference implementation.
    #mpow = mvec .^ degree
    mpow = copy(mvec)
    for i in 1:degree-1
        mpow .*= mvec
    end

    cipher = encrypt(rng, enc_key, mvec, log_precision, log_cap)

    cipher_res = power_of_2(mul_key, cipher, log_precision, log_degree)

    dvec = decrypt(secret_key, cipher_res)

    test_approx(mpow, dvec, 20)
end


@testcase "power()" begin
    n = 2^6
    degree = 13
    log_precision = 30
    log_cap = 300

    rng = MyRNG(12345)
    params = Params(log_polynomial_length=8, log_lo_modulus=300)

    secret_key = SecretKey(rng, params)

    enc_key = EncryptionKey(rng, secret_key)
    mul_key = MultiplicationKey(rng, secret_key)

    mvec = [randomCircle(rng) for i in 1:n] # randn(rng, n) + im * randn(rng, n)

    # Temporary replacement for bit-to-bit compatibility with the reference implementation.
    #mpow = mvec .^ degree
    mpow = copy(mvec)
    for i in 1:degree-1
        mpow .*= mvec
    end

    cipher = encrypt(rng, enc_key, mvec, log_precision, log_cap)

    cipher_res = power(mul_key, cipher, log_precision, degree)

    dvec = decrypt(secret_key, cipher_res)

    test_approx(mpow, dvec, 20)
end


@testcase "inv()" begin
    n = 2^6
    log_precision = 30
    log_cap = 300
    steps = 6

    rng = MyRNG(12345)
    params = Params(log_polynomial_length=8, log_lo_modulus=300)

    secret_key = SecretKey(rng, params)

    enc_key = EncryptionKey(rng, secret_key)
    mul_key = MultiplicationKey(rng, secret_key)

    # TODO: if one uses 1.0 instead of 0.1, the error becomes huge. Why?
    mvec = [randomCircle(rng, 0.1) for i in 1:n] # randn(rng, n) + im * randn(rng, n)
    minv = 1 ./ mvec

    cipher = encrypt(rng, enc_key, mvec, log_precision, log_cap)

    cipher_res = inv(mul_key, cipher, log_precision, steps)

    dvec = decrypt(secret_key, cipher_res)

    test_approx(minv, dvec, 22)
end


@testcase "log_plus_one()" begin
    n = 2^6
    log_precision = 30
    log_cap = 300
    degree = 7

    rng = MyRNG(12345)
    params = Params(log_polynomial_length=8, log_lo_modulus=300)

    secret_key = SecretKey(rng, params)

    enc_key = EncryptionKey(rng, secret_key)
    mul_key = MultiplicationKey(rng, secret_key)

    mvec = [randomComplex(rng, 0.1) for i in 1:n] # randn(rng, n) + im * randn(rng, n)
    mlog = log.(mvec .+ 1)

    cipher = encrypt(rng, enc_key, mvec, log_precision, log_cap)

    cipher_res = log_plus_one(mul_key, cipher, log_precision, degree)

    dvec = decrypt(secret_key, cipher_res)

    test_approx(mlog, dvec, 23)
end


@testcase "exp()" begin
    n = 2^6
    log_precision = 30
    log_cap = 300
    degree = 7

    rng = MyRNG(12345)
    params = Params(log_polynomial_length=8, log_lo_modulus=300)

    secret_key = SecretKey(rng, params)

    enc_key = EncryptionKey(rng, secret_key)
    mul_key = MultiplicationKey(rng, secret_key)

    mvec = [randomComplex(rng) for i in 1:n] # randn(rng, n) + im * randn(rng, n)
    mexp = exp.(mvec)

    cipher = encrypt(rng, enc_key, mvec, log_precision, log_cap)

    cipher_res = exp(mul_key, cipher, log_precision, degree)

    dvec = decrypt(secret_key, cipher_res)

    test_approx(mexp, dvec, 12)
end


@testcase "div_by_po2()" begin
    n = 2^6
    log_precision = 30
    log_cap = 300
    bits = 20

    rng = MyRNG(12345)
    params = Params(log_polynomial_length=8, log_lo_modulus=300)

    secret_key = SecretKey(rng, params)

    enc_key = EncryptionKey(rng, secret_key)

    mvec = [randomComplex(rng) for i in 1:n] # randn(rng, n) + im * randn(rng, n)
    mdiv = mvec ./ 2^bits

    cipher = encrypt(rng, enc_key, mvec, log_precision, log_cap)

    cipher_res = div_by_po2(cipher, bits)

    dvec = decrypt(secret_key, cipher_res)

    test_approx(mdiv, dvec, 23)
end


@testcase "sigmoid()" for lazy in ["lazy", "eager"]
    n = 2^6
    log_precision = 30
    log_cap = 300
    degree = 7

    rng = MyRNG(12345)
    params = Params(log_polynomial_length=8, log_lo_modulus=300)

    secret_key = SecretKey(rng, params)

    enc_key = EncryptionKey(rng, secret_key)
    mul_key = MultiplicationKey(rng, secret_key)

    mvec = [randomComplex(rng) for i in 1:n] # randn(rng, n) + im * randn(rng, n)
    msig = exp.(mvec) ./ (1 .+ exp.(mvec))

    cipher = encrypt(rng, enc_key, mvec, log_precision, log_cap)

    cipher_res = sigmoid(mul_key, cipher, log_precision, degree, lazy=(lazy == "lazy"))

    dvec = decrypt(secret_key, cipher_res)

    test_approx(msig, dvec, 12)
end


@testcase "bootstrap" begin
    rng = MyRNG(12345)

    log_t = 4 # means that we use Taylor approximation in [-1/t,1/t] with double angle fomula
    log_n = 3
    n = 2^log_n
    log_precision = 20
    log_cap = 30
    # TODO (see issue #1): why are we setting boot context precision as `log_cap + 4`?
    bc_precision = log_cap + 4

    params = Params(log_polynomial_length=8, log_lo_modulus=600)
    log_plen = params.log_polynomial_length

    secret_key = SecretKey(rng, params)
    enc_key = EncryptionKey(rng, secret_key)
    mul_key = MultiplicationKey(rng, secret_key)
    conj_key = ConjugationKey(rng, secret_key)

    bk = BootstrapKey(rng, secret_key, enc_key, mul_key, conj_key, log_n, bc_precision)

    mvec = randomComplexArray(rng, n)

    cipher = encrypt(rng, bk.enc_key, mvec, log_precision, log_cap)

    cipher_res = bootstrap(bk, cipher, log_t)

    dvec = decrypt(secret_key, cipher_res)

    # The main purpose of bootstrapping
    @test cipher_res.log_cap - cipher_res.log_precision > cipher.log_cap - cipher.log_precision

    test_approx(mvec, dvec, 6)
end


end
