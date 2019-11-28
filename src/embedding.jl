struct EmbeddingPlan
    max_len :: Int
    rot_group :: Array{Int, 1}
    ksi_pows :: Array{Complex{Float64}, 1}

    function EmbeddingPlan(max_len::Int)
        rot_group = Array{Int}(undef, max_len ÷ 2)
        five_pows = 1
        for i in 1:max_len÷2
            rot_group[i] = five_pows
            five_pows *= 5
            five_pows = mod(five_pows, max_len * 2)
        end

        ksi_pows = Array{Complex{Float64}}(undef, max_len * 2 + 1)
        m_pi = Float64(pi)
        for j in 1:max_len*2
            angle = 2.0 * m_pi * (j-1) / (max_len * 2)
            ksi_pows[j] = mycos(angle) + im * mysin(angle) # exp(im * angle)
        end
        ksi_pows[max_len * 2 + 1] = ksi_pows[1]

        new(max_len, rot_group, ksi_pows)
    end
end


function arrayBitReverse!(vals::Array{Complex{Float64}, 1}, n::Int)
    j = 0
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


function unembed(plan::EmbeddingPlan, vals::Array{Complex{Float64}, 1})

    # TODO: assuming n is a power of 2
    n = length(vals)
    log_n = trailing_zeros(n)
    new_vals = copy(vals)

    for log_len in log_n:-1:0
        len = 2^log_len
        for i in 0:len:n-1
            lenh = len >> 1
            lenq = len << 2
            gap = plan.max_len * 2 ÷ lenq
            for j in 0:lenh-1
                idx = (lenq - (plan.rot_group[j+1] % lenq)) * gap
                u = new_vals[i + j + 1] + new_vals[i + j + lenh + 1]
                v = new_vals[i + j + 1] - new_vals[i + j + lenh + 1]
                v *= plan.ksi_pows[idx + 1]
                new_vals[i + j + 1] = u
                new_vals[i + j + lenh + 1] = v
            end
        end
    end

    arrayBitReverse!(new_vals, n)

    for i in 0:n-1
        new_vals[i+1] /= n
    end
    new_vals
end


function embed(plan::EmbeddingPlan, vals::Array{Complex{Float64}, 1})
    new_vals = copy(vals)
    n = length(vals)
    arrayBitReverse!(new_vals, n)

    # TODO: assuming n is a power of 2
    log_n = trailing_zeros(n)

    for log_len in 1:log_n
        len = 1 << log_len

        for i in 0:len:n-1
            lenh = len >> 1
            lenq = len << 2
            gap = plan.max_len * 2 ÷ lenq
            for j in 0:lenh-1
                idx = (plan.rot_group[j+1] % lenq) * gap
                u = new_vals[i + j + 1]
                v = new_vals[i + j + lenh + 1]
                v *= plan.ksi_pows[idx + 1]
                new_vals[i + j + 1] = u + v
                new_vals[i + j + lenh + 1] = u - v
            end
        end
    end

    new_vals
end


# TODO: replace with a macro?
const _embedding_plans = IdDict{Params, EmbeddingPlan}()


function embedding_plan(params::Params)
    if haskey(_embedding_plans, params)
        _embedding_plans[params]
    else
        res = EmbeddingPlan(2^params.log_polynomial_length)
        _embedding_plans[params] = res
        res
    end
end


function embed(params::Params, vals::Array{Complex{Float64}, 1})
    embed(embedding_plan(params), vals)
end


function unembed(params::Params, vals::Array{Complex{Float64}, 1})
    unembed(embedding_plan(params), vals)
end

