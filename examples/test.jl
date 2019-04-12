push!(LOAD_PATH, "../src")

using HEAAN
using HEAAN:
    Ring, SecretKey, MyRNG, Scheme, randomComplex, encryptSingle, decryptSingle,
    randomComplexArray, encrypt, hexfloat, decrypt,
    add, mult, imult,
    addLeftRotKey!, leftRotateFast,
    addConjKey!, conjugate,
    randomCircle, SchemeAlgo, powerOf2,
    power,
    inverse


function testEncryptSingle(logq::Int, logp::Int)

    rng = MyRNG(12345)

    ring = Ring()
    secretKey = SecretKey(rng, ring)
    scheme = Scheme(rng, secretKey, ring)

    mval = randomComplex(rng)
    cipher = encryptSingle(rng, scheme, mval, logp, logq)

    dval = decryptSingle(scheme, secretKey, cipher)

    println("mval=$mval")
    println("dval=$dval")
    println("diff=$(mval-dval)")
end


function testEncrypt(logq::Int, logp::Int, logn::Int)

    rng = MyRNG(12345)

    ring = Ring()
    secretKey = SecretKey(rng, ring)
    scheme = Scheme(rng, secretKey, ring)

    n = 2^logn
    mvec = randomComplexArray(rng, n)
    cipher = encrypt(rng, scheme, mvec, n, logp, logq)

    dvec = decrypt(scheme, secretKey, cipher)

    for i in 1:n
        println(hexfloat(real(mvec[i])), " ", hexfloat(real(dvec[i])))
        println(hexfloat(imag(mvec[i])), " ", hexfloat(imag(dvec[i])))
    end
end


function testAdd(logq::Int, logp::Int, logn::Int)

    rng = MyRNG(12345)
    ring = Ring()
    secretKey = SecretKey(rng, ring)
    scheme = Scheme(rng, secretKey, ring)

    n = 2^logn
    mvec1 = randomComplexArray(rng, n)
    mvec2 = randomComplexArray(rng, n)
    madd = mvec1 .+ mvec2

    cipher1 = encrypt(rng, scheme, mvec1, n, logp, logq)
    cipher2 = encrypt(rng, scheme, mvec2, n, logp, logq)

    cipher_res = add(scheme, cipher1, cipher2)

    dadd = decrypt(scheme, secretKey, cipher_res)

    for i in 1:n
        println(hexfloat(real(madd[i])), " ", hexfloat(real(dadd[i])))
        println(hexfloat(imag(madd[i])), " ", hexfloat(imag(dadd[i])))
    end
end


function testMult(logq::Int, logp::Int, logn::Int)

    rng = MyRNG(12345)
    ring = Ring()
    secretKey = SecretKey(rng, ring)
    scheme = Scheme(rng, secretKey, ring)

    n = 2^logn
    mvec1 = randomComplexArray(rng, n)
    mvec2 = randomComplexArray(rng, n)
    mmul = mvec1 .* mvec2

    cipher1 = encrypt(rng, scheme, mvec1, n, logp, logq)
    cipher2 = encrypt(rng, scheme, mvec2, n, logp, logq)

    cipher_res = mult(scheme, cipher1, cipher2)

    dmul = decrypt(scheme, secretKey, cipher_res)

    for i in 1:n
        println(hexfloat(real(mmul[i])), " ", hexfloat(real(dmul[i])))
        println(hexfloat(imag(mmul[i])), " ", hexfloat(imag(dmul[i])))
    end
end


function testimult(logq::Int, logp::Int, logn::Int)

    rng = MyRNG(12345)
    ring = Ring()
    secretKey = SecretKey(rng, ring)
    scheme = Scheme(rng, secretKey, ring)

    n = 2^logn
    mvec = randomComplexArray(rng, n)
    imvec = mvec .* im

    cipher = encrypt(rng, scheme, mvec, n, logp, logq)
    cipher_res = imult(scheme, cipher)

    idvec = decrypt(scheme, secretKey, cipher_res)

    for i in 1:n
        println(hexfloat(real(imvec[i])), " ", hexfloat(real(idvec[i])))
        println(hexfloat(imag(imvec[i])), " ", hexfloat(imag(idvec[i])))
    end
end


function testRotateFast(logq::Int, logp::Int, logn::Int, logr::Int)

    rng = MyRNG(12345)

    ring = Ring()
    secretKey = SecretKey(rng, ring)
    scheme = Scheme(rng, secretKey, ring)

    n = 2^logn
    r = 2^logr

    addLeftRotKey!(rng, scheme, secretKey, r)

    mvec = randomComplexArray(rng, n)
    cipher = encrypt(rng, scheme, mvec, n, logp, logq)

    cipher_res = leftRotateFast(scheme, cipher, r)

    dvec = decrypt(scheme, secretKey, cipher_res)

    mvec_rot = circshift(mvec, -r)

    for i in 1:n
        println(hexfloat(real(mvec_rot[i])), " ", hexfloat(real(dvec[i])))
        println(hexfloat(imag(mvec_rot[i])), " ", hexfloat(imag(dvec[i])))
    end
end


function testConjugate(logq::Int, logp::Int, logn::Int)

    rng = MyRNG(12345)

    ring = Ring()
    secretKey = SecretKey(rng, ring)
    scheme = Scheme(rng, secretKey, ring)

    n = 2^logn

    addConjKey!(rng, scheme, secretKey)

    mvec = randomComplexArray(rng, n)
    mvecconj = conj.(mvec)

    cipher = encrypt(rng, scheme, mvec, n, logp, logq)

    cipher_res = conjugate(scheme, cipher)

    dvecconj = decrypt(scheme, secretKey, cipher_res)

    for i in 1:n
        println(hexfloat(real(mvecconj[i])), " ", hexfloat(real(dvecconj[i])))
        println(hexfloat(imag(mvecconj[i])), " ", hexfloat(imag(dvecconj[i])))
    end
end


function testPowerOf2(logq::Int, logp::Int, logn::Int, logdeg::Int)

    rng = MyRNG(12345)

    ring = Ring()
    secretKey = SecretKey(rng, ring)
    scheme = Scheme(rng, secretKey, ring)
    algo = SchemeAlgo(scheme)

    n = 2^logn
    degree = 2^logdeg

    mvec = [randomCircle(rng) for i in 1:n]
    # Temporary replacement for bit-to-bit compatibility with the reference implementation.
    #mpow = mvec .^ degree
    mpow = copy(mvec)
    for i in 1:degree-1
        mpow .*= mvec
    end

    cipher = encrypt(rng, scheme, mvec, n, logp, logq)

    cpow = powerOf2(algo, cipher, logp, logdeg)

    dpow = decrypt(scheme, secretKey, cpow)

    for i in 1:n
        println(hexfloat(real(mpow[i])), " ", hexfloat(real(dpow[i])))
        println(hexfloat(imag(mpow[i])), " ", hexfloat(imag(dpow[i])))
    end
end



function testPower(logq::Int, logp::Int, logn::Int, degree::Int)

    rng = MyRNG(12345)

    ring = Ring()
    secretKey = SecretKey(rng, ring)
    scheme = Scheme(rng, secretKey, ring)
    algo = SchemeAlgo(scheme)

    n = 2^logn

    mvec = [randomCircle(rng) for i in 1:n]
    # Temporary replacement for bit-to-bit compatibility with the reference implementation.
    #mpow = mvec .^ degree
    mpow = copy(mvec)
    for i in 1:degree-1
        mpow .*= mvec
    end

    cipher = encrypt(rng, scheme, mvec, n, logp, logq)

    cpow = power(algo, cipher, logp, degree)

    dpow = decrypt(scheme, secretKey, cpow)

    for i in 1:n
        println(hexfloat(real(mpow[i])), " ", hexfloat(real(dpow[i])))
        println(hexfloat(imag(mpow[i])), " ", hexfloat(imag(dpow[i])))
    end
end


function testInverse(logq::Int, logp::Int, logn::Int, steps::Int)

    rng = MyRNG(12345)

    ring = Ring()
    secretKey = SecretKey(rng, ring)
    scheme = Scheme(rng, secretKey, ring)
    algo = SchemeAlgo(scheme)

    n = 2^logn

    mvec = [randomCircle(rng, 0.1) for i in 1:n]
    minv = 1 ./ mvec

    cipher = encrypt(rng, scheme, mvec, n, logp, logq)

    cinv = inverse(algo, cipher, logp, steps)

    dinv = decrypt(scheme, secretKey, cinv)

    for i in 1:n
        println(hexfloat(real(minv[i])), " ", hexfloat(real(dinv[i])))
        println(hexfloat(imag(minv[i])), " ", hexfloat(imag(dinv[i])))
    end
end

#testEncryptSingle(100, 30)
#testEncrypt(100, 30, 4)
#testAdd(100, 30, 4)
#testMult(100, 30, 4)
#testimult(100, 30, 4)
#testRotateFast(100, 30, 4, 1)
#testConjugate(100, 30, 4)
#testPowerOf2(300, 30, 4, 4)
#testPower(300, 30, 4, 13)
testInverse(300, 25, 4, 6)
