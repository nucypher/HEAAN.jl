module HEAAN

using Random
using Primes
using DarkIntegers


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


mutable struct MyRNG
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


include("params.jl")
include("ntl.jl")
include("evaluator_utils.jl")
include("ring_multiplier.jl")
include("bootcontext.jl")
include("ring.jl")
include("secret_key.jl")
include("key.jl")
include("plaintext.jl")
include("ciphertext.jl")
include("scheme.jl")
include("scheme_algo.jl")

end
