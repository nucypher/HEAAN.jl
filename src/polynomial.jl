struct RNSPolynomial
    plan :: RNSPlan
    residuals :: Array{UInt64, 2}
end


struct RNSPolynomialTransformed
    plan :: RNSPlan
    residuals :: Array{UInt64, 2}
end


function _to_rns(plan::RNSPlan, x::Polynomial{BinModuloInt{T, Q}}, np::Int) where {T, Q}
    plen = length(x.coeffs)
    res = Array{UInt64}(undef, plen, np)
    for i in 1:plen
        res[i,:] .= to_rns_signed(plan, x.coeffs[i], np)
    end
    RNSPolynomial(plan, res)
end


function _ntt_rns(plan::RNSPlan, x::Array{UInt64, 2}, inverse::Bool)
    np = size(x, 2)
    res = similar(x)
    for j in 1:np
        r = x[:,j]
        m = plan.primes[j]

        # TODO: keep residuals already casted to RRElem?
        tp = RRElem{UInt64, m}
        rr = DarkIntegers.ntt(tp.(r), inverse=inverse, negacyclic=true)

        res[:,j] .= DarkIntegers.rr_value.(rr)
    end
    res
end


function _ntt_forward(x::RNSPolynomial)
    RNSPolynomialTransformed(x.plan, _ntt_rns(x.plan, x.residuals, false))
end


function _ntt_inverse(x::RNSPolynomialTransformed)
    RNSPolynomial(x.plan, _ntt_rns(x.plan, x.residuals, true))
end


function to_rns_transformed(plan::RNSPlan, x::Polynomial{BinModuloInt{T, Q}}, np::Int) where {T, Q}
    _ntt_forward(_to_rns(plan, x, np))
end

to_rns_transformed(plan::RNSPlan, x::Polynomial) = to_rns_transformed(plan, x, length(plan.primes))


function _from_rns(::Type{Polynomial{BinModuloInt{T, Q}}}, x::RNSPolynomial) where {T, Q}
    plan = x.plan
    plen = size(x.residuals, 1)
    res = Array{BinModuloInt{T, Q}}(undef, plen)
    for i in 1:plen
        res[i] = from_rns_signed(plan, BinModuloInt{T, Q}, x.residuals[i,:])
    end
    Polynomial(res, true)
end


function from_rns_transformed(x::RNSPolynomialTransformed, log_modulus::Int)
    tp = Polynomial{BinModuloInt{BigInt, log_modulus}}
    _from_rns(tp, _ntt_inverse(x))
end


function Base.:+(x::RNSPolynomialTransformed, y::RNSPolynomialTransformed)
    # TODO: pick the minimum number of residuals for `x` and `y` for the result?
    plan = x.plan
    res = similar(x.residuals)
    np = size(x.residuals, 2)
    for j in 1:np
        res[:,j] = addmod.(x.residuals[:,j], y.residuals[:,j], plan.primes[j])
    end
    RNSPolynomialTransformed(plan, res)
end


function Base.:*(x::RNSPolynomialTransformed, y::RNSPolynomialTransformed)
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
    RNSPolynomialTransformed(plan, res)
end


function mult(x::Polynomial{BinModuloInt{T, Q}}, y::RNSPolynomialTransformed, np::Int) where {T, Q}
    plan = y.plan
    x_rns = to_rns_transformed(plan, x, np)
    from_rns_transformed(x_rns * y, Q)
end


# TODO: seems to correspond to Rotate() in the paper?
# Or does it refer to `circshift` itself?
# "For an input encryption of m(Y), return an encryption of m(Y^(5^k)) in the same level"
function left_rotate(x::Polynomial, r::Integer)
    res = Polynomial(similar(x.coeffs), x.negacyclic)
    n = length(x.coeffs)
    # This corresponds to the rotation group used in embedding.jl
    # Essentially it's rotation_group[r], but this function is too low-level to
    # create EmbeddingPlan here.
    pow = mod(5^r, 2n)
    for i in 0:n-1
        shift = mod(i * pow, 2n)
        if shift < n
            res.coeffs[shift+1] = x.coeffs[i+1]
        else
            res.coeffs[shift - n + 1] = -x.coeffs[i+1]
        end
    end
    res
end


function conjugate(x::Polynomial)
    res = Polynomial(similar(x.coeffs), x.negacyclic)
    res.coeffs[1] = x.coeffs[1]
    res.coeffs[2:end] .= .-x.coeffs[end:-1:2]
    res
end


function mod_down_by(x::Polynomial{BinModuloInt{T, Q}}, dq::Int) where {T, Q}
    Polynomial(trunc.(BinModuloInt{T, Q - dq}, x.coeffs), x.negacyclic)
end


function mod_down_to(x::Polynomial{BinModuloInt{T, Q}}, log_q::Int) where {T, Q}
    Polynomial(trunc.(BinModuloInt{T, log_q}, x.coeffs), x.negacyclic)
end


function mod_up_to(x::Polynomial{BinModuloInt{T, Q}}, log_q::Int) where {T, Q}
    # TODO: used only in bootstrap
    # Two ways are possible: just increase the range,
    # or keep negative numbers negative (that is x -> x - Q_old + Q_new).
    # Not sure what is the right one.
    Polynomial(convert.(BinModuloInt{T, log_q}, x.coeffs), x.negacyclic)
end


# TODO: do we even need this as a separate function?
function mul_by_monomial(x::Polynomial, pwr::Integer)
    shift_polynomial(x, pwr)
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
