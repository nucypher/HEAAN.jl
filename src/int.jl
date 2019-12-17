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


struct BinModuloInt{T, Q} <: Unsigned
    value :: T

    function BinModuloInt{T, Q}(x::T) where {T, Q}
        new{T, Q}(x)
    end

    function BinModuloInt(x::T, log_modulus::Int) where T
        # TODO: (issue #11) check that log_modulus actually fits in T
        new{T, log_modulus}(x)
    end
end


@inline Base.promote_type(::Type{BinModuloInt{T, Q}}, ::Type{BinModuloInt{T, Q}}) where {T, Q} =
    BinModuloInt{T, Q}
@inline Base.promote_type(::Type{BinModuloInt{T, Q}}, ::Type{<:Integer}) where {T, Q} =
    BinModuloInt{T, Q}
@inline Base.promote_type(::Type{<:Integer}, ::Type{BinModuloInt{T, Q}}) where {T, Q} =
    BinModuloInt{T, Q}


function num_bits(x::BinModuloInt)
    if signbit(x)
        x = -x
    end
    num_bits(x.value)
end


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
    BinModuloInt{T, Q}(normalize(convert(T, x), Q))


function mod_down_to(::Type{BinModuloInt{T, Q1}}, x::BinModuloInt{T, Q2}) where {T, Q1, Q2}
    @assert Q1 <= Q2
    BinModuloInt{T, Q1}(x.value & all_bits_mask(BinModuloInt{T, Q1}))
end


function Base.convert(::Type{BinModuloInt{T, Q1}}, x::BinModuloInt{T, Q2}) where {T, Q1, Q2}
    @assert Q1 >= Q2
    neg = signbit(x)
    BinModuloInt{T, Q1}(neg ? x.value - modulus(T, Q2) + modulus(T, Q1) : x.value)
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


function Base.:>>(x::BinModuloInt{T, Q}, shift::Integer) where {T, Q}
    BinModuloInt{T, Q - shift}(x.value >> shift)
end


#=
What C++ HEAAN is doing is essentially

    right_shift_rounded(x, shift) = (signbit(x) ? -1 : 1) * (abs(x + 2^(shift-1)) >> shift)

This doesn't quite correspond to `round(x / 2^shift)`, which is the intention,
and reduces precision. The Julia version of that would be

    half_shift = BinModuloInt{T, Q}(one(T) << (shift - 1))
    new_x = x + half_shift
    (signbit(x) && signbit(new_x)) ? -((-new_x) >> shift) : new_x >> shift

Instead, we're implementing an algorithm that does correspond to `round(x / 2^shift)`
(except for when the result is `-2^(Q-shift-1)`).
=#
function right_shift_rounded(x::BinModuloInt{T, Q}, shift::Integer) where {T, Q}

    s = signbit(x)
    x = s ? -x : x
    v = x.value + modulus(T, shift - 1)
    tie = iszero(v & all_bits_mask(T, shift)) # that is, whether x/2^shift = ###.5
    v >>= shift

    # Using the "Round half to even" tie-breaking rule, the same as Julia uses by default.
    if tie && isodd(v)
        v -= one(T)
    end

    # Our integer range is `(-2^(q-1), 2^(q-1)]`, where `q = Q - shift`.
    # So we have to deal with the value `-2^(q-1)` somehow.
    # One way is to wrap it around to `2^(q-1)`, but that will mean that we get a positive
    # number by rounding a negative number, which can lead to big errors.
    # Instead, we just add 1, which leads to very small and rare errors.
    if s && v == one(T) << (Q - shift - 1)
        v -= one(T)
    end

    v = HEAAN.BinModuloInt{T, Q - shift}(v)
    s ? -v : v
end


Base.:<<(x::BinModuloInt{T, Q}, shift::Int) where {T, Q} =
    BinModuloInt{T, Q + shift}(x.value << shift)


Base.string(x::BinModuloInt) = x.value
Base.show(io::IO, x::BinModuloInt) = print(io, string(x))
