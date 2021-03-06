"""
    SecretKey(rng::AbstractRNG, params::Params)

Secret key, used for decryption.
Takes a [`Params`](@ref) object.
"""
struct SecretKey

    params :: Params
    nonzero_entries :: Array{Pair{Int, Bool}, 1} # a list of (power, is_minus_one)

    function SecretKey(rng::AbstractRNG, params::Params)
        # Sample from HWT distribution in the paper:
        # {1, -1} at `secret_key_length` positions, the rest are zeros.
        positions = sample_unique(
            rng, 0, 2^params.log_polynomial_length-1, params.secret_key_length)
        bits = rand(rng, Bool, params.secret_key_length)

        # We are recording the key in sparse form
        new(params, [pos => bit for (pos, bit) in zip(positions, bits)])
    end
end


function sample_unique(rng::AbstractRNG, a::Int, b::Int, len::Int)
    res = Set{Int}()
    while length(res) < len
        push!(res, rand(rng, a:b))
    end
    sort(collect(res))
end


function Base.:*(secret_key::SecretKey, pa::Polynomial{BinModuloInt{T, Q}}) where {T, Q}

    sk_list = secret_key.nonzero_entries

    power, minus_one = sk_list[1]
    res = mul_by_monomial(pa, power)
    if minus_one
        res = -res
    end

    for (power, minus_one) in sk_list[2:end]
        if minus_one
            res -= mul_by_monomial(pa, power)
        else
            res += mul_by_monomial(pa, power)
        end
    end

    res
end


function as_polynomial(secret_key::SecretKey, log_modulus::Int)
    params = secret_key.params
    tp = BinModuloInt{BigInt, log_modulus}
    sk_poly = Polynomial(zeros(tp, 2^params.log_polynomial_length), negacyclic_modulus)
    for (power, minus_one) in secret_key.nonzero_entries
        sk_poly.coeffs[power+1] = minus_one ? -one(tp) : one(tp)
    end
    sk_poly
end


function square(secret_key::SecretKey, log_modulus::Int)
    secret_key * as_polynomial(secret_key, log_modulus)
end
