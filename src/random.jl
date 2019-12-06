function rand_big_int(rng::AbstractRNG, log_modulus::Int, dims...)
    coeffs = rand(rng, zero(BigInt):high_bit_mask(BigInt, log_modulus), dims...)
    Polynomial(BinModuloInt{BigInt, log_modulus}.(coeffs), true)
end


function sample_ZO(rng::AbstractRNG, len::Int, log_modulus::Int)
    zero_val = zero(BinModuloInt{BigInt, log_modulus})
    one_val = one(BinModuloInt{BigInt, log_modulus})
    minus_one = -one_val

    bits1 = rand(rng, Bool, len)
    bits2 = rand(rng, Bool, len)
    res = [
        b1 ? zero_val : (b2 ? one_val : minus_one)
        for (b1, b2) in zip(bits1, bits2)]

    Polynomial(res, true)
end


function discrete_gaussian(rng::AbstractRNG, stddev::Float64, dims...)
    round.(Int, randn(rng, dims...) .* stddev)
end