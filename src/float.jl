significand_bits(::Type{T}) where T <: Base.IEEEFloat = precision(T) - 1
significand_bits(::T) where T <: Base.IEEEFloat = significand_bits(T)


# There's already a built-in `significand()` that does not return what we want,
# so calling this one `mantissa()`.
mantissa(x::T) where T <: Base.IEEEFloat =
    reinterpret(Unsigned, x) & Base.significand_mask(T)


exponent_bias(::Type{T}) where T <: Base.IEEEFloat =
    Int(Base.exponent_one(T) >> significand_bits(T))


"""
Returns `x * 2^shift` rounded to an integer.
If `x` is negative, return `2^log_modulus - abs(x) * 2^shift`.
The resulting value is guaranteed to lie in `[0, 2^log_modulus)` range,
which means that `-(2^shift-1)/(2^(log_modulus-1)) <= x <= 2^shift/(2^(log_modulus-1))`.

Negative values are encoded in the high half of the interval `[0, 2^log_modulus)`
(the sign of `x` will correspond to what `is_negative()` returns for the resulting integer).
"""
function float_to_integer(
        ::Type{V}, x::T, shift::Int, log_modulus::Int) where {V <: Integer, T <: Base.IEEEFloat}

    if iszero(x)
        # `exponent()` does not work for x=0
        return zero(V)
    end

    m = mantissa(x)
    sb = significand_bits(x)
    u_type = typeof(m)

    # Add the implied 1 to the fractional part
    m |= one(u_type) << sb

    # Shift for the mantissa bits.
    # Since mantissa is the fractional part, its length should be subtracted
    mantissa_shift = exponent(x) + shift - sb

    # At this point x == m * 2.0^shift

    correction = if mantissa_shift < 0
        # We will have a fractional part after the shift, need to round

        # The fractional part (that will be discarded after the shift),
        # that is the lowest (-shift) bits of `m`
        remainder = m & ((one(u_type) << (-mantissa_shift)) - one(u_type))

        # See if it's >= 0.5
        leading_zeros(remainder) == (sizeof(u_type) * 8 - (-mantissa_shift))
    else
        false
    end

    # Convert the mantissa and shift it
    r = convert(V, m) << mantissa_shift

    if correction
        r += one(V)
    end

    n = num_bits(r)
    s = signbit(x)

    # Check the range. A special check for the maximum positive integer.
    if n >= log_modulus && !(!s && r == high_bit_mask(V, log_modulus))
        throw(DomainError(
            x,
            "the value must lie in range " *
            "[-(2^shift-1)/(2^(log_modulus-1)), 2^shift/(2^(log_modulus-1))]"))
    end

    s ? modulus(V, log_modulus) - r : r
end


function float_to_integer(::Type{V}, x::BigFloat, shift::Int, log_modulus::Int) where V <: Integer
    xc = copy(x)
    xc.exp += shift
    xi = round(V, xc)
    res = convert(V, xi)
    signbit(xi) ? modulus(V, log_modulus) + res : res
end


"""
Convert an integer `x` to float and divide by `2^shift`.
`x` must lie in range `[0, 2^log_modulus)`.

Negative values are encoded in the high half of the interval `[0, 2^log_modulus)`
(`is_negative()` is used to determine the sign).
"""
function integer_to_float(
        ::Type{V}, x::T, shift::Int, log_modulus::Int) where {T <: Integer, V <: Base.IEEEFloat}

    if iszero(x)
        return zero(V)
    end

    n = num_bits(x)
    if signbit(x) || n > log_modulus
        throw(DomainError(x, "the value must lie in range [0, 2^log_modulus)"))
    end

    s = is_negative(x, log_modulus)
    # TODO: move to a function?
    if s
        x = modulus(T, log_modulus) - x
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


float_to_integer(::Type{BinModuloInt{T, Q}}, x, shift::Int) where {T, Q} =
    BinModuloInt{T, Q}(float_to_integer(T, x, shift, Q))
integer_to_float(tp, x::BinModuloInt{T, Q}, shift::Int) where {T, Q} =
    integer_to_float(tp, x.value, shift, Q)
