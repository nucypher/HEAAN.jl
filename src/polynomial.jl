struct RNSPolynomial
    plan :: RNSPlan
    residuals :: Array{UInt64, 2}
end


function to_rns(plan::RNSPlan, x::Polynomial{BinModuloInt{T, Q}}, np::Int) where {T, Q}
    plen = length(x.coeffs)
    res = Array{UInt64}(undef, plen, np)
    for i in 1:plen
        res[i,:] .= to_rns_signed(plan, x.coeffs[i], np)
    end
    RNSPolynomial(plan, res)
end

# TODO: `length(plan.primes)` should be the default?
to_rns(plan::RNSPlan, x::Polynomial) = to_rns(plan, x, length(plan.primes))


function from_rns(::Type{Polynomial{BinModuloInt{T, Q}}}, x::RNSPolynomial, np::Int) where {T, Q}
    plan = x.plan
    plen = size(x.residuals, 1)
    res = Array{BinModuloInt{T, Q}}(undef, plen)
    for i in 1:plen
        res[i] = from_rns_signed(plan, BinModuloInt{T, Q}, x.residuals[i,:])
    end
    Polynomial(res, true)
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


function Base.:+(x::RNSPolynomial, y::RNSPolynomial)
    # TODO: pick the minimum number of residuals for `x` and `y` for the result?
    plan = x.plan
    res = similar(x.residuals)
    np = size(x.residuals, 2)
    for j in 1:np
        res[:,j] = addmod.(x.residuals[:,j], y.residuals[:,j], plan.primes[j])
    end
    RNSPolynomial(plan, res)
end


function Base.:*(x::RNSPolynomial, y::RNSPolynomial)
    # TODO: keep keys and such in M-representation, to speed up multiplication
    # (although make sure `+` is still processed correctly)
    # TODO: pick the minimum number of residuals for `x` and `y` for the result?
    plan = x.plan
    res = similar(x.residuals)
    np = size(x.residuals, 2)
    for j in 1:np
        # TODO: use Barrett reduction, or Montgomery multiplication
        res[:,j] = mulmod.(x.residuals[:,j], y.residuals[:,j], plan.primes[j])
    end
    RNSPolynomial(plan, res)
end


function mult(x::T, y::RNSPolynomial, np::Int) where T <: Polynomial
    plan = y.plan
    x_rns = ntt_rns(to_rns(plan, x, np), inverse=false)
    res_rns = ntt_rns(x_rns * y, inverse=true)
    from_rns(T, res_rns, np)
end


# FIXME: these methods will not be necessary when Polynomial broadcasting
# is implemented in DarkIntegers


Base.:+(x::AbstractArray, y::Polynomial{BinModuloInt{T, Q}}) where {T, Q} =
    Polynomial(convert.(BinModuloInt{T, Q}, x), y.negacyclic) + y


Base.:-(x::AbstractArray, y::Polynomial{BinModuloInt{T, Q}}) where {T, Q} =
    Polynomial(convert.(BinModuloInt{T, Q}, x), y.negacyclic) - y


Base.:>>(x::Polynomial{BinModuloInt{T, Q}}, shift::Integer) where {T, Q} =
    Polynomial(x.coeffs .>> shift, x.negacyclic)

Base.:<<(x::Polynomial{BinModuloInt{T, Q}}, shift::Integer) where {T, Q} =
    Polynomial(x.coeffs .<< shift, x.negacyclic)
