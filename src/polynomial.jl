struct RNSPolynomial
    plan :: RNSPlan
    residuals :: Array{UInt64, 2}

    # The numbers stored are guaranteed to fit into `log_modulus` bits,
    # that is they lie in `(P - 2^(log_modulus-1) + 1, P - 1] U [0, 2^(log_modulus-1)]`,
    # where `P = prod(plan.primes[1:np])`, and `np` is the second dimension of `residuals`
    # (and, of course, `P > 2^log_modulus - 1`, so these intervals don't intersect).
    log_modulus :: Int

    # The maximum possible range for stored numbers for the given number of residuals.
    # `log_modulus` can't get greater than this.
    # `log_range` has a one-to-one correspondence with the number of primes used.
    log_range :: Int

    polynomial_modulus :: DarkIntegers.AbstractPolynomialModulus

    function RNSPolynomial(plan, residuals, log_modulus, log_range, polynomial_modulus)
        @assert log_modulus <= log_range
        @assert log_range == max_log_modulus(plan, size(residuals, 2))
        new(plan, residuals, log_modulus, log_range, polynomial_modulus)
    end
end


struct RNSPolynomialTransformed
    plan :: RNSPlan
    residuals :: Array{UInt64, 2}
    log_modulus :: Int
    log_range :: Int
    polynomial_modulus :: DarkIntegers.AbstractPolynomialModulus

    function RNSPolynomialTransformed(plan, residuals, log_modulus, log_range, polynomial_modulus)
        @assert log_modulus <= log_range
        new(plan, residuals, log_modulus, log_range, polynomial_modulus)
    end
end


function _to_rns(plan::RNSPlan, x::Polynomial{BinModuloInt{T, Q}}, log_range::Int) where {T, Q}

    np = min_nprimes(plan, log_range)
    log_range = max_log_modulus(plan, np)

    plen = length(x.coeffs)
    res = Array{UInt64}(undef, plen, np)
    for i in 1:plen
        res[i,:] .= to_rns_signed(plan, x.coeffs[i], np)
    end
    RNSPolynomial(plan, res, Q, log_range, x.modulus)
end


function _ntt_rns(plan::RNSPlan, x::Array{UInt64, 2}, inverse::Bool, negacyclic::Bool)
    np = size(x, 2)
    res = similar(x)
    for j in 1:np
        r = x[:,j]
        m = plan.primes[j]

        # TODO: (issue #9) keep residuals already casted to ModUInt?
        tp = ModUInt{UInt64, m}
        rr = DarkIntegers.ntt(tp.(r), inverse=inverse, tangent=negacyclic)

        res[:,j] .= value.(rr)
    end
    res
end


function _ntt_forward(x::RNSPolynomial)
    RNSPolynomialTransformed(
        x.plan,
        _ntt_rns(x.plan, x.residuals, false, x.polynomial_modulus == negacyclic_modulus),
        x.log_modulus,
        x.log_range,
        x.polynomial_modulus)
end


function _ntt_inverse(x::RNSPolynomialTransformed)
    RNSPolynomial(
        x.plan,
        _ntt_rns(x.plan, x.residuals, true, x.polynomial_modulus == negacyclic_modulus),
        x.log_modulus,
        x.log_range,
        x.polynomial_modulus)
end


function to_rns_transformed(
        plan::RNSPlan, x::Polynomial{BinModuloInt{T, Q}}, add_range::Int=0) where {T, Q}

    # The intention is to have enough range for one multiplication of this polynomial
    # and another one with `log_modulus <= add_range`.
    log_range = _log_modulus_mul(Q, add_range, length(x.coeffs))

    _ntt_forward(_to_rns(plan, x, log_range))
end

to_rns_transformed(plan::RNSPlan, x::Polynomial) = to_rns_transformed(plan, x, length(plan.primes))


function _from_rns(::Type{Polynomial{BinModuloInt{T, Q}}}, x::RNSPolynomial) where {T, Q}
    @assert Q <= x.log_range
    plan = x.plan
    plen = size(x.residuals, 1)
    res = Array{BinModuloInt{T, Q}}(undef, plen)
    for i in 1:plen
        res[i] = from_rns_signed(plan, BinModuloInt{T, Q}, x.residuals[i,:])
    end
    Polynomial(res, x.polynomial_modulus)
end


function from_rns_transformed(x::RNSPolynomialTransformed, log_modulus::Int=0)
    if log_modulus == 0
        log_modulus = x.log_modulus
    else
        @assert log_modulus <= x.log_modulus
    end
    tp = Polynomial{BinModuloInt{BigInt, log_modulus}}
    _from_rns(tp, _ntt_inverse(x))
end


function Base.:+(x::RNSPolynomialTransformed, y::RNSPolynomialTransformed)

    @assert x.polynomial_modulus == y.polynomial_modulus

    # The modulus can grow at worst by one bit
    # (if we're adding two maximum or two minimum possible numbers)
    new_log_modulus = max(x.log_modulus, y.log_modulus) + 1

    new_log_range = min(x.log_range, y.log_range)
    @assert new_log_modulus <= new_log_range

    plan = x.plan
    np = min(size(x.residuals, 2), size(y.residuals, 2))
    res = Array{UInt64}(undef, size(x.residuals, 1), np)
    for j in 1:np
        res[:,j] = addmod.(x.residuals[:,j], y.residuals[:,j], plan.primes[j])
    end
    RNSPolynomialTransformed(plan, res, new_log_modulus, new_log_range, x.polynomial_modulus)
end


function _log_modulus_mul(log_modulus1::Int, log_modulus2::Int, polynomial_length::Int)
    # The result is the product of polynomials, so the maximum number we can encounter
    # is the polynomial-length sum of products of maximum numbers from `x` and `y`
    # (perhaps that's even a bit too conservative for negacyclic polynomials).
    log_plen = num_bits(polynomial_length) - 1

    # `-1` is because each coefficient <= 2^(log_modulus-1)
    # So the maximum is <= 2^(log_modulus1-1) * 2^(log_modulus2-1) * 2^log_plen
    # `+1` because we need to fit both positive and negative numbers of that range.
    (log_modulus1 - 1) + (log_modulus2 - 1) + log_plen + 1
end


function Base.:*(x::RNSPolynomialTransformed, y::RNSPolynomialTransformed)
    # TODO: (issue #9) keep keys and such in M-representation, to speed up multiplication
    # (although make sure `+` is still processed correctly)

    @assert x.polynomial_modulus == y.polynomial_modulus

    plen = size(x.residuals, 1)
    new_log_modulus = _log_modulus_mul(x.log_modulus, y.log_modulus, plen)

    new_log_range = min(x.log_range, y.log_range)
    @assert new_log_modulus <= new_log_range

    plan = x.plan
    np = min(size(x.residuals, 2), size(y.residuals, 2))
    res = Array{UInt64}(undef, size(x.residuals, 1), np)
    for j in 1:np
        # TODO: (issue #9) use Barrett reduction, or Montgomery multiplication
        res[:,j] = mulmod.(x.residuals[:,j], y.residuals[:,j], plan.primes[j])
    end
    RNSPolynomialTransformed(plan, res, new_log_modulus, new_log_range, x.polynomial_modulus)
end


poly_lm(x::Polynomial{BinModuloInt{T, Q}}) where {T, Q} = Q


function mul_by_rns(x::Polynomial{BinModuloInt{T, Q}}, y::RNSPolynomialTransformed) where {T, Q}
    @assert x.modulus == y.polynomial_modulus
    plan = y.plan
    x_rns = _ntt_forward(_to_rns(plan, x, y.log_range))
    from_rns_transformed(x_rns * y, Q)
end


function left_rotate(x::Polynomial, r::Integer)
    res = Polynomial(similar(x.coeffs), x.modulus)
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
    res = Polynomial(similar(x.coeffs), x.modulus)
    res.coeffs[1] = x.coeffs[1]
    res.coeffs[2:end] .= .-x.coeffs[end:-1:2]
    res
end


function mod_down_by(x::Polynomial{BinModuloInt{T, Q}}, dq::Int) where {T, Q}
    Polynomial(mod_down_to.(BinModuloInt{T, Q - dq}, x.coeffs), x.modulus)
end


function mod_down_to(x::Polynomial{BinModuloInt{T, Q}}, log_q::Int) where {T, Q}
    Polynomial(mod_down_to.(BinModuloInt{T, log_q}, x.coeffs), x.modulus)
end


function mod_up_to(x::Polynomial{BinModuloInt{T, Q}}, log_q::Int) where {T, Q}
    Polynomial(convert.(BinModuloInt{T, log_q}, x.coeffs), x.modulus)
end
