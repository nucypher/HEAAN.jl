@generated function _all_bits_mask(::Type{T}, ::Val{Q}) where {T, Q}
    mask = (one(T) << Q) - one(T)
    :( $mask )
end


"""
Returns the precomputed value `(one(T) << log_modulus) - one(T)`.
"""
all_bits_mask(::Type{T}, log_modulus::Int) where T = _all_bits_mask(T, Val(log_modulus))


@generated function _modulus(::Type{T}, ::Val{Q}) where {T, Q}
    modulus = one(T) << Q
    :( $modulus )
end


"""
Returns the precomputed value `one(T) << log_modulus`.
"""
modulus(::Type{T}, log_modulus::Int) where T = _modulus(T, Val(log_modulus))


@generated function _high_bit_mask(::Type{T}, ::Val{Q}) where {T, Q}
    mask = one(T) << (Q - 1)
    :( $mask )
end


"""
Returns the precomputed value `one(T) << (log_modulus - 1)`.
"""
high_bit_mask(::Type{T}, log_modulus::Int) where T = _high_bit_mask(T, Val(log_modulus))


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


struct BinModuloInt{T, Q} <: Number
    value :: T

    function BinModuloInt{T, Q}(x::T) where {T, Q}
        new{T, Q}(x)
    end

    function BinModuloInt(x::T, log_modulus::Int) where T
        # TODO: check that log_modulus actually fits in T
        new{T, log_modulus}(x)
    end
end


num_bits(x::BinModuloInt) = num_bits(x.value)


all_bits_mask(::Type{BinModuloInt{T, Q}}) where {T, Q} = all_bits_mask(T, Q)
modulus(::Type{BinModuloInt{T, Q}}) where {T, Q} = modulus(T, Q)
high_bit_mask(::Type{BinModuloInt{T, Q}}) where {T, Q} = high_bit_mask(T, Q)
Base.signbit(x::BinModuloInt{T, Q}) where {T, Q} = is_negative(x.value, Q)


function normalize(x::T, log_modulus::Int) where T
    x & all_bits_mask(T, log_modulus)
end


Base.zero(::Type{BinModuloInt{T, Q}}) where {T, Q} = BinModuloInt{T, Q}(zero(T))


Base.one(::Type{BinModuloInt{T, Q}}) where {T, Q} = BinModuloInt{T, Q}(one(T))


Base.convert(::Type{BinModuloInt{T, Q}}, x::Integer) where {T, Q} =
    BinModuloInt{T, Q}(normalize(big(x), Q))


# TODO: check that it works correctly with negative numbers
# TODO: actually, `trunc` is not the right function to use here - for numbers is assumes that
# the number is in the final range already.
# Need our own function.
function Base.trunc(::Type{BinModuloInt{T, Q1}}, x::BinModuloInt{T, Q2}) where {T, Q1, Q2}
    @assert Q1 <= Q2
    BinModuloInt{T, Q1}(x.value & all_bits_mask(BinModuloInt{T, Q1}))
end

function Base.convert(::Type{BinModuloInt{T, Q1}}, x::BinModuloInt{T, Q2}) where {T, Q1, Q2}
    @assert Q1 >= Q2
    BinModuloInt{T, Q1}(x.value)
end


Base.:+(x::BinModuloInt{T, Q}, y::BinModuloInt{T, Q}) where {T, Q} =
    BinModuloInt(normalize(x.value + y.value, Q), Q)


Base.:-(x::BinModuloInt{T, Q}, y::BinModuloInt{T, Q}) where {T, Q} =
    BinModuloInt(normalize(x.value - y.value, Q), Q)

Base.:-(x::BinModuloInt{T, Q}) where {T, Q} =
    BinModuloInt(iszero(x.value) ? x.value : modulus(BinModuloInt{T, Q}) - x.value, Q)


Base.:*(x::BinModuloInt{T, Q}, y::BinModuloInt{T, Q}) where {T, Q} =
    BinModuloInt(normalize(x.value * y.value, Q), Q)


Base.:(==)(x::BinModuloInt{T, Q}, y::BinModuloInt{T, Q}) where {T, Q} =
    x.value == y.value


# TODO: check that it works correctly with all corner cases
Base.:>>(x::BinModuloInt{T, Q}, shift::Integer) where {T, Q} =
    BinModuloInt{T, Q - shift}((x.value + (one(T) << (shift - 1))) >> shift)

Base.:<<(x::BinModuloInt{T, Q}, shift::Integer) where {T, Q} =
    BinModuloInt{T, Q + shift}(x.value << shift)


Base.string(x::BinModuloInt) = x.value
Base.show(io::IO, x::BinModuloInt) = print(io, string(x))
