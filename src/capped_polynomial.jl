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
    x & all_bits_mask(T, log_modulus)
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
        q = modulus(BigInt, log_modulus)
        q - ((q - x) >> shift)
    else
        x >> shift
    end
end

# TODO: check if rounding/not rounding makes a difference for precision retaining
Base.:>>(x::CappedPolynomial{T, Q}, shift::Integer) where {T, Q} =
    CappedPolynomial{T, Q - shift}(Polynomial((x.polynomial.coeffs .+ (one(BigInt) << (shift-1))) .>>  shift, true))



struct RNSPolynomial
    plan :: RNSPlan
    residuals :: Array{UInt64, 2}
end


function to_rns(plan::RNSPlan, x::CappedPolynomial{T, Q}, np::Int) where {T, Q}
    plen = length(x.polynomial.coeffs)
    res = Array{UInt64}(undef, plen, np)
    for i in 1:plen
        res[i,:] .= to_rns(plan, x.polynomial.coeffs[i], np, Q)
    end
    RNSPolynomial(plan, res)
end

# TODO: `length(plan.primes)` should be the default?
to_rns(plan::RNSPlan, x::CappedPolynomial{T, Q}) where {T, Q} = to_rns(plan, x, length(plan.primes))


function from_rns(::Type{CappedPolynomial{T, Q}}, x::RNSPolynomial, np::Int) where {T, Q}
    plan = x.plan
    plen = size(x.residuals, 1)
    res = Array{BigInt}(undef, plen)
    for i in 1:plen
        res[i] = from_rns(plan, x.residuals[i,:], Q)
    end
    CappedPolynomial{T, Q}(Polynomial(res, true))
end


function ntt_rns(x::RNSPolynomial; inverse::Bool=false)
    plan = x.plan
    res = Array{UInt64}(undef, size(x.residuals)...)
    np = size(x.residuals, 2)
    for j in 1:np
        r = x.residuals[:,j]
        m = plan.primes[j]

        # TODO: keep residuals already casted to RRElem?
        tp = RRElem{UInt64, m}
        rr = DarkIntegers.ntt(tp.(r), inverse=inverse, negacyclic=true)

        res[:,j] .= DarkIntegers.rr_value.(rr)
    end

    RNSPolynomial(plan, res)
end


function Base.:*(x::RNSPolynomial, y::RNSPolynomial)
    # TODO: keep keys and such in M-representation, to speed up multiplication
    plan = x.plan
    res = similar(x.residuals)
    np = size(x.residuals, 2)
    for j in 1:np
        # TODO: use Barrett reduction, or Montgomery multiplication
        res[:,j] = mulmod.(x.residuals[:,j], y.residuals[:,j], plan.primes[j])
    end
    RNSPolynomial(plan, res)
end


function mult(x::CappedPolynomial{T, Q}, y::RNSPolynomial, np::Int) where {T, Q}
    plan = y.plan
    x_rns = ntt_rns(to_rns(plan, x, np), inverse=false)
    res_rns = ntt_rns(x_rns * y, inverse=true)
    from_rns(CappedPolynomial{T, Q}, res_rns, np)
end
