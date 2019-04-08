push!(LOAD_PATH, "../src")

using HEAAN
using HEAAN: Ring, SecretKey, MyRNG, Scheme, randomComplex, encryptSingle, decryptSingle


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


testEncryptSingle(100, 30)
