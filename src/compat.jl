mysin(x) = ccall((:sin, "libc"), Float64, (Float64,), x)
mycos(x) = ccall((:cos, "libc"), Float64, (Float64,), x)


function hexfloat(x::Float64)
    ux = reinterpret(UInt64, x)
    m = ux & ((1 << 52) - 1)
    ux = (ux - m) >> 52
    e = signed(ux & ((1 << 11) - 1)) - 1023
    s = ux & 0x800

    ms = string(m, base=16)
    if length(ms) < 13
        ms = "0"^(13 - length(ms)) * ms
    end
    ms = replace(ms, r"0+$" => s"")

    ss = s == 0 ? "" : "-"
    es = e < 0 ? "" : "+"
    "$(ss)0x1.$(ms)p$(es)$e"
end


mutable struct MyRNG <: AbstractRNG
    state :: UInt64
end


function myrand(rng::MyRNG)
    x = rng.state # The state must be seeded with a nonzero value
    x = xor(x, (x >> 12)) # a
    x = xor(x, (x << 25)) # b
    x = xor(x, (x >> 27)) # c
    rng.state = x
    x * UInt64(0x2545F4914F6CDD1D)
end


function myrand_float(rng::MyRNG)
    r = myrand(rng)
    convert(Float64, r) / 18446744073709551616.
end


function myRandomBits_ZZ(rng::MyRNG, len::Int)
    res = zero(BigInt)
    shift = 0
    while true
        u64 = myrand(rng)
        bitlen = len > 64 ? 64 : len
        for i in 0:bitlen-1
            if (u64 & (one(UInt64) << i)) != 0
                res += one(BigInt) << (shift + i)
            end
        end
        shift += 64
        len -= 64
        if len < 0
            break
        end
    end
    res
end



# Intended replacement for `myRandomBits_ZZ`, since what we need is essentially
# just an array of random bits.
function myRandomBits_array(rng::MyRNG, len::Int)
    res = [false for i in 1:len]
    shift = 0
    while true
        u64 = myrand(rng)
        bitlen = len > 64 ? 64 : len
        for i in 0:bitlen-1
            if (u64 & (one(UInt64) << i)) != 0
                res[shift + i + 1] = true
            end
        end
        shift += 64
        len -= 64
        if len < 0
            break
        end
    end
    res
end



function myRandomBits_long(rng::MyRNG, len::Int)
    res = zero(Int)
    shift = 0
    while true
        u64 = myrand(rng)
        bitlen = len > 64 ? 64 : len
        for i in 0:bitlen-1
            if (u64 & (one(UInt64) << i)) != 0
                res += one(Int) << (shift + i)
            end
        end
        shift += 64
        len -= 64
        if len < 0
            break
        end
    end
    res
end


function myrand_long(rng::MyRNG, lim::Int)
    Int(mod(myrand(rng), lim))
end


# TODO: should be a part of DarkIntegers, and use `mulmod()`
function powmod(x::T, y::Integer, m::T) where T
    @assert y >= 0
    if y == 0
        one(T)
    elseif y == 1
        x
    else
        acc = one(T)
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


randomReal(rng::MyRNG, bound::Float64 = 1.0) = myrand_float(rng) * bound
randomComplex(rng::MyRNG, bound::Float64 = 1.0) =
    randomReal(rng, bound) + im * randomReal(rng, bound)
randomComplexArray(rng::MyRNG, n::Int, bound::Float64 = 1.0) =
    [randomComplex(rng, bound) for i in 1:n]


function randomCircle(rng::MyRNG, anglebound::Float64 = 1.0)
    angle_ = randomReal(rng, anglebound)
    M_PI = Float64(pi)
    mycos(angle_ * 2 * M_PI) + im * mysin(angle_ * 2 * M_PI)
end


function sampleUniform2(rng::MyRNG, bits::Int, N::Int)
    [myRandomBits_ZZ(rng, bits) for i in 1:N]
end


function rand_gauss(rng::MyRNG, N::Int, sigma::Float64)

    # x is in range [0, q)

    bignum = Int(0xfffffff)

    res = Array{Float64}(undef, N)

    for i in 0:2:N-1
        # TODO: essentially Box-Muller sampling
        r1 = (1 + myrand_long(rng, bignum)) / (bignum + 1)
        r2 = (1 + myrand_long(rng, bignum)) / (bignum + 1)
        theta = 2pi * r1
        rr = sqrt(-2.0 * log(r2)) * sigma

        res[i+1] = rr * mycos(theta)
        res[i+2] = rr * mysin(theta)
    end

    res
end
