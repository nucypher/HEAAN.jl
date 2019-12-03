struct SecretKey

    params :: Params
    nonzero_entries :: Array{Pair{Int, Bool}, 1} # a list of (position, is_minus_one)

    function SecretKey(rng::AbstractRNG, params::Params)
        # Sample from HWT distribution in the paper:
        # {1, -1} at `secret_key_length` positions, the rest are zeros.

        #positions = sample_unique(rng, 1, 2^params.log_polynomial_length, params.secret_key_length)
        #bits = rand(rng, Bool, params.secret_key_length)

        # We are recording the key in sparse form
        #new(params, collect(zip(positions, bits)))

        res = Dict{Int, Bool}()
        idx = 0
        tmp = myRandomBits_ZZ(rng, params.secret_key_length)
        while idx < params.secret_key_length
            i = myRandomBits_long(rng, params.log_polynomial_length)
            if !haskey(res, i+1)
                res[i+1] = !(tmp & (one(BigInt) << idx) == 0)
                idx += 1
            end
        end

        new(params, collect(res))
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
    res = shift_polynomial(pa, power - 1)
    if minus_one
        res = -res
    end

    for (power, minus_one) in sk_list[2:end]
        if minus_one
            res -= shift_polynomial(pa, power - 1)
        else
            res += shift_polynomial(pa, power - 1)
        end
    end

    res
end


function as_polynomial(secret_key::SecretKey)
    params = secret_key.params
    tp = BinModuloInt{BigInt, params.log_lo_modulus}
    sk_poly = Polynomial(zeros(tp, 2^params.log_polynomial_length), true)
    for (power, minus_one) in secret_key.nonzero_entries
        sk_poly.coeffs[power] = minus_one ? -one(tp) : one(tp)
    end
    sk_poly
end


function square(secret_key::SecretKey)
    as_polynomial(secret_key) * secret_key
end
