significand_bits(::Type{T}) where T <: Base.IEEEFloat = precision(T) - 1
significand_bits(::T) where T <: Base.IEEEFloat = significand_bits(T)


# There's already a built-in `significand()` that does not return what we want,
# so calling this one `mantissa()`.
mantissa(x::T) where T <: Base.IEEEFloat =
    reinterpret(Unsigned, x) & Base.significand_mask(T)


"""
Returns `x * 2^shift` rounded to an integer.
"""
function float_to_integer(::Type{V}, x::T, shift::Int, log_full::Int) where {V <: Integer, T <: Base.IEEEFloat}

    if iszero(x)
        # `exponent()` does not work for x=0
        return zero(V)
    end

    m = mantissa(x)
    sb = significand_bits(x)
    u_type = typeof(m)

    # Add the implied 1 to the fractional part
    m |= one(u_type) << sb

    # Total shift.
    # Since mantissa is the fractional part, its length should be subtracted
    full_shift = exponent(x) + shift - sb

    # At this point x == m * 2.0^shift

    correction = if full_shift < 0
        # We will have a fractional part after the shift, need to round

        # The fractional part (that will be discarded after the shift),
        # that is the lowest (-shift) bits of `m`
        remainder = m & ((one(u_type) << (-full_shift)) - one(u_type))

        # See if it's >= 0.5
        leading_zeros(remainder) == (sizeof(u_type) * 8 - (-full_shift))
    else
        false
    end

    # Convert the mantissa and shift it
    r = V(m) << full_shift

    if correction
        r += one(V)
    end

    @assert num_bits(r) < log_full

    signbit(x) ? (one(BigInt) << log_full) - r : r
end


function float_to_integer(::Type{V}, x::BigFloat, shift::Int, log_full::Int) where V <: Integer
    xc = copy(x)
    xc.exp += shift
    xi = round(BigInt, xc)
    res = convert(V, xi)
    signbit(xi) ? (one(BigInt) << log_full) + res : res
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


exponent_bias(::Type{T}) where T <: Base.IEEEFloat =
    Int(Base.exponent_one(T) >> significand_bits(T))


"""
Convert an integer `x` to float and divide by `2^shift`.
"""
function integer_to_float(::Type{V}, x::T, shift::Int, log_modulus::Int) where {T <: Integer, V <: Base.IEEEFloat}

    if iszero(x)
        return zero(V)
    end

    s = is_negative(x, log_modulus)
    # TODO: move to a function?
    if s
        x = (one(T) << log_modulus) - x
    end

    n = num_bits(x)
    e = n - 1 - shift + exponent_bias(V)
    if e < 0
        return zero(V)
    end

    s_len = significand_bits(V)
    u_type = Base.uinttype(V)

    # `temp` now contains one high bit which we are going to drop (since it's implied in a float),
    # `s_len` mantissa bits and one low bit which tells us whether we need to round up or down.
    temp = convert(u_type, x >> (n - s_len - 2))

    # Getting rid of the low bit and the high bit
    s_bits = (temp >> 1) & Base.significand_mask(V)

    # FIXME: For some reason the rounding in the reference function
    # (Float64(BigFloat(x) / 2^shift)) works in a slightly strange way.
    # For remainders strictly less or greater than `b.1000...000` it works as expected:
    # rounds down or up, respectively.
    # But for the `remainder == b.1000...000`, it rounds up if the last bit of the mantissa is 1
    # and down otherwise.
    #
    # So in order to preserve compatibility we have to check for that.
    # Can probably be removed later.

    first_remainder_bit = isodd(temp)
    last_mantissa_bit = isodd(s_bits)
    remainder_trailing_zeros = trailing_zeros(x)
    middle_remainder = remainder_trailing_zeros == n - s_len - 2

    if first_remainder_bit && (!middle_remainder || last_mantissa_bit)
        s_bits += 1
    end

    e_bits = convert(u_type, e) << s_len
    res = reinterpret(V, s_bits | e_bits)

    s ? -res : res
end
