module HEAAN

using Random
using Primes
using DarkIntegers


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
include("ring.jl")
include("secret_key.jl")
include("key.jl")
include("plaintext.jl")
include("ciphertext.jl")
include("scheme.jl")

end
