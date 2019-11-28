
@generated function all_bits_mask(::Type{T}, ::Val{Q}) where {T, Q}
    mask = (one(T) << Q) - one(T)
    :( $mask )
end

@generated function modulus(::Type{T}, ::Val{Q}) where {T, Q}
    modulus = one(T) << Q
    :( $modulus )
end

#=
@generated function high_bit_mask(::Type{BinModuloInt{T, Q}}) where {T, Q}
    mask = one(T) << (Q - 1)
    :( $mask )
end


=#


function is_negative(x::BigInt, logq::Int)
    # In the paper the range of numbers is (-q/2, q/2]
    # So we consider everything in [0, q/2] positive and [q/2+1, q-1) negative.
    half_q = one(BigInt) << (logq - 1)
    high_bit = !iszero(x & half_q)
    high_bit && x != half_q
end


struct CappedPolynomial{T <: Polynomial, Q}
    polynomial :: T

    function CappedPolynomial{T, Q}(polynomial::T) where {T, Q}
        #check_range(coeffs, logq)
        new{T, Q}(polynomial)
    end

    function CappedPolynomial(polynomial::T, log_modulus::Int) where T
        #check_range(coeffs, logq)
        new{T, log_modulus}(polynomial)
    end
end


function log_modulus(::CappedPolynomial{T, Q}) where {T, Q}
    Q
end


function normalize(x::T, log_modulus::Int) where T
    x & all_bits_mask(T, Val(log_modulus))
end


function normalize(x::Polynomial{T}, log_modulus::Int) where T
    # TODO: return CappedPolynomial here? Or just make it one of the constructors?
    Polynomial(normalize.(x.coeffs, log_modulus), true)
end


Base.:+(x::CappedPolynomial{T, Q}, y::CappedPolynomial{T, Q}) where {T, Q} =
    CappedPolynomial{T, Q}(normalize(x.polynomial + y.polynomial, Q))

Base.:+(x::Array{Int, 1}, y::CappedPolynomial{T, Q}) where {T, Q} =
    CappedPolynomial{T, Q}(normalize(Polynomial(convert.(BigInt, x), true), Q)) + y

Base.:-(x::CappedPolynomial{T, Q}, y::CappedPolynomial{T, Q}) where {T, Q} =
    CappedPolynomial{T, Q}(normalize(x.polynomial - y.polynomial, Q))

Base.:-(x::Array{Int, 1}, y::CappedPolynomial{T, Q}) where {T, Q} =
    CappedPolynomial{T, Q}(normalize(Polynomial(convert.(BigInt, x), true), Q)) - y

Base.:*(x::CappedPolynomial{T, Q}, y::CappedPolynomial{T, Q}) where {T, Q} =
    CappedPolynomial{T, Q}(normalize(x.polynomial * y.polynomial), Q)


function my_rshift(x::BigInt, shift::Integer, log_modulus::Int)
    if is_negative(x, log_modulus)
        q = modulus(BigInt, Val(log_modulus))
        q - ((q - x) >> shift)
    else
        x >> shift
    end
end

# TODO: check if rounding/not rounding makes a difference for precision retaining
Base.:>>(x::CappedPolynomial{T, Q}, shift::Integer) where {T, Q} =
    CappedPolynomial{T, Q - shift}(Polynomial((x.polynomial.coeffs .+ (one(BigInt) << (shift-1))) .>>  shift, true))



struct RNSPolynomial
    rns :: RNS
    residuals :: Array{UInt64, 2}
end


function to_rns(rns::RNS, x::CappedPolynomial{T, Q}, np::Int) where {T, Q}
    plen = length(x.polynomial.coeffs)
    res = Array{UInt64}(undef, plen, np)
    for i in 1:plen
        res[i,:] .= to_rns(rns, x.polynomial.coeffs[i], np, Q)
    end
    RNSPolynomial(rns, res)
end


function from_rns(::Type{CappedPolynomial{T, Q}}, x::RNSPolynomial, np::Int) where {T, Q}
    rns = x.rns
    plen = size(x.residuals, 1)
    res = Array{BigInt}(undef, plen)
    for i in 1:plen
        res[i] = from_rns(rns, x.residuals[i,:], np, Q)
    end
    CappedPolynomial{T, Q}(Polynomial(res, true))
end


function ntt_rns(x::RNSPolynomial; inverse::Bool=false)
    rns = x.rns
    res = Array{UInt64}(undef, size(x.residuals)...)
    np = size(x.residuals, 2)
    for j in 1:np
        r = x.residuals[:,j]
        m = rns.pVec[j]

        # TODO: keep residuals already casted to RRElem?
        tp = RRElem{UInt64, m}
        rr = DarkIntegers.ntt(tp.(r), inverse=inverse, negacyclic=true)

        res[:,j] .= DarkIntegers.rr_value.(rr)
    end

    RNSPolynomial(rns, res)
end


function Base.:*(x::RNSPolynomial, y::RNSPolynomial)
    # TODO: keep keys and such in M-representation, to speed up multiplication
    rns = x.rns
    res = similar(x.residuals)
    np = size(x.residuals, 2)
    for j in 1:np
        res[:,j] = mulmod.(x.residuals[:,j], y.residuals[:,j], rns.pVec[j])
    end
    RNSPolynomial(rns, res)
end


function mult(x::CappedPolynomial{T, Q}, y::RNSPolynomial, np::Int) where {T, Q}
    rns = y.rns
    x_rns = ntt_rns(to_rns(rns, x, np), inverse=false)
    res_rns = ntt_rns(x_rns * y, inverse=true)
    from_rns(CappedPolynomial{T, Q}, res_rns, np)
end
