@generated function all_bits_mask(::Type{T}, ::Val{Q}) where {T, Q}
    mask = (one(T) << Q) - one(T)
    :( $mask )
end


"""
Returns the precomputed value `(one(T) << log_modulus) - one(T)`.
"""
all_bits_mask(::Type{T}, log_modulus::Int) where T = all_bits_mask(T, Val(log_modulus))


@generated function modulus(::Type{T}, ::Val{Q}) where {T, Q}
    modulus = one(T) << Q
    :( $modulus )
end


"""
Returns the precomputed value `one(T) << log_modulus`.
"""
modulus(::Type{T}, log_modulus::Int) where T = modulus(T, Val(log_modulus))


@generated function high_bit_mask(::Type{T}, ::Val{Q}) where {T, Q}
    mask = one(T) << (Q - 1)
    :( $mask )
end


"""
Returns the precomputed value `one(T) << (log_modulus - 1)`.
"""
high_bit_mask(::Type{T}, log_modulus::Int) where T = high_bit_mask(T, Val(log_modulus))


"""
For `x` in range `[0, 2^log_modulus)`, returns `true` if `x > 2^log_modulus`
and `false` otherwise.
"""
function is_negative(x::T, log_modulus::Int) where T
    # In the paper the range of numbers is (-modulus/2, modulus/2]
    # So we consider everything in [0, modulus/2] positive and [modulus/2+1, modulus-1) negative.
    high_bit = !iszero(x >> (log_modulus - 1))
    high_bit && x != high_bit_mask(T, log_modulus)
end


"""
Return the number of bits in the representation of the absolute value of an integer.
"""
num_bits(x::T) where T <: Unsigned = DarkIntegers.bitsizeof(T) - leading_zeros(x)

function num_bits(x::T) where T <: Signed
    if x == typemin(T)
        throw(DomainError(
            x,
            "typemin() for fixed-size signed integers " *
            "is not supported in (-q/2,q/2] representation"))
    end

    num_bits(unsigned(signbit(x) ? -x : x))
end

function num_bits(x::BigInt)
    if iszero(x)
        return 0
    end

    if signbit(x)
        x = -x
    end

    for i in abs(x.size):-1:1
        limb = unsafe_load(x.d, i)

        # BigInts seem to resize automatically, but we're venturing
        # into the undocumented territory here, so just in case
        # handle possible empty limbs.
        if !iszero(limb)
            limb_bits = sizeof(limb) << 3
            return (i - 1) * limb_bits + (limb_bits - leading_zeros(limb))
        end
    end
end
