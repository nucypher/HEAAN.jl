struct EmbeddingPlan
    max_len :: Int
    rotation_group :: Array{Int, 1}
    root_powers :: Array{Complex{Float64}, 1}

    function EmbeddingPlan(max_len::Int)

        rotation_group = mod.(5 .^ (0:max_len-1), max_len * 4)

        # The last element being the same as the first one is intentional -
        # it simplifies addressing slightly.
        root_powers = exp.(2pi * im .* (0:max_len*4) / (max_len * 4))

        new(max_len, rotation_group, root_powers)
    end
end


"""
Permute the given array as `x[i] -> x[bitreverse(i)]`, `i = 0...n-1`.

`n` must be a power of 2.
"""
function array_bit_reverse!(vals::Array{Complex{Float64}, 1})
    j = 0
    n = length(vals)
    for i in 1:n-1
        bit = n >> 1
        while j >= bit
            j -= bit
            bit >>= 1
        end
        j += bit
        if i < j
            t = vals[i+1]
            vals[i+1] = vals[j+1]
            vals[j+1] = t
        end
    end
end


"""
Calculate the canonical embedding.
Essentially, for a polynomial of `n-1`-th power `a(x)` (given as an array of coefficients),
calculate `n` values `[a(xi), a(xi^5), a(xi^25), ...]` where `xi = exp(2pi * im / (4n))`.

`n` must be a power of 2.
"""
function embed(plan::EmbeddingPlan, vals::Array{Complex{Float64}, 1})

    n = length(vals)
    log_n = trailing_zeros(n)

    new_vals = copy(vals)

    array_bit_reverse!(new_vals)

    for log_len in 1:log_n
        len = 1 << log_len

        for i in 0:len:n-1
            lenh = len >> 1
            lenq = len << 2
            gap = plan.max_len * 4 ÷ lenq
            for j in 0:lenh-1
                idx = (plan.rotation_group[j+1] % lenq) * gap
                u = new_vals[i + j + 1]
                v = new_vals[i + j + lenh + 1]
                v *= plan.root_powers[idx + 1]
                new_vals[i + j + 1] = u + v
                new_vals[i + j + lenh + 1] = u - v
            end
        end
    end

    new_vals
end


"""
An operation inverse to `embed()`. For the given array of evaluations of the polynomial
finds the polynomial itself.
"""
function unembed(plan::EmbeddingPlan, vals::Array{Complex{Float64}, 1})

    n = length(vals)
    log_n = trailing_zeros(n)

    new_vals = copy(vals)

    for log_len in log_n:-1:0
        len = 1 << log_len
        for i in 0:len:n-1
            lenh = len >> 1
            lenq = len << 2
            gap = plan.max_len * 4 ÷ lenq
            for j in 0:lenh-1
                idx = (lenq - (plan.rotation_group[j+1] % lenq)) * gap
                u = new_vals[i + j + 1] + new_vals[i + j + lenh + 1]
                v = new_vals[i + j + 1] - new_vals[i + j + lenh + 1]
                v *= plan.root_powers[idx + 1]
                new_vals[i + j + 1] = u
                new_vals[i + j + lenh + 1] = v
            end
        end
    end

    array_bit_reverse!(new_vals)

    new_vals ./= n

    new_vals
end
