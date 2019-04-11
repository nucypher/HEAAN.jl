push!(LOAD_PATH, "../src")

using HEAAN
using HEAAN:
    Ring, SecretKey, MyRNG, Scheme, randomComplex, encryptSingle, decryptSingle,
    randomComplexArray, encrypt, hexfloat, decrypt


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


#testEncryptSingle(100, 30)
testEncrypt(100, 30, 4)
