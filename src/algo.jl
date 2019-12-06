function power_bin_exp(mk::MultiplicationKey, cipher::Ciphertext, log_precision::Int, log_pwr::Int)
    res = cipher
    for i in 1:log_pwr
        res = square(mk, res)
        res = rescale_by(res, log_precision)
    end
    res
end


function power(mk::MultiplicationKey, cipher::Ciphertext, log_precision::Int, pwr::Int)
    log_pwr = num_bits(pwr) - 1
    res = power_bin_exp(mk, cipher, log_precision, log_pwr)
    rem_pwr = pwr - (1 << log_pwr)
    if rem_pwr > 0
        tmp = power(mk, cipher, log_precision, rem_pwr)
        bits_down = tmp.log_cap - res.log_cap
        tmp = mod_down_by(tmp, bits_down)
        res = mul(mk, res, tmp)
        res = rescale_by(res, log_precision)
    end
    res
end


function Base.inv(mk::MultiplicationKey, cipher::Ciphertext, log_precision::Int, steps::Int)
    cbar = negate(cipher)
    cbar = add_const(cbar, 1.0)
    cpow = cbar
    tmp = add_const(cbar, 1.0)
    tmp = mod_down_by(tmp, log_precision)
    res = tmp
    for i in 1:steps-1
        cpow = square(mk, cpow)
        cpow = rescale_by(cpow, log_precision)
        tmp = add_const(cpow, 1.0)
        tmp = mul(mk, tmp, res)
        tmp = rescale_by(tmp, log_precision)
        res = tmp
    end
    res
end


function power_of_2_extended(
        mk::MultiplicationKey, cipher::Ciphertext, log_precision::Int, log_degree::Int)
    res = Array{Ciphertext}(undef, log_degree + 1)
    res[1] = cipher
    for i in 1:log_degree
        res[i+1] = rescale_by(square(mk, res[i]), log_precision)
    end
    res
end


function power_extended(mk::MultiplicationKey, cipher::Ciphertext, log_precision::Int, degree::Int)
    log_degree = floor(Int, log2(degree))
    cpows = power_of_2_extended(mk, cipher, log_precision, log_degree)
    res = Array{Ciphertext}(undef, degree)
    idx = 0
    for i in 0:log_degree-1
        powi = 1 << i
        res[idx+1] = cpows[i+1]
        idx += 1
        for j in 0:powi-1-1
            bits_down = res[j+1].log_cap - cpows[i+1].log_cap
            res[idx+1] = mod_down_by(res[j+1], bits_down)
            res[idx+1] = mul(mk, res[idx+1], cpows[i+1])
            res[idx+1] = rescale_by(res[idx+1], log_precision)
            idx += 1
        end
    end
    res[idx+1] = cpows[log_degree+1]
    idx += 1
    degree2 = 1 << log_degree
    for i in 0:degree-degree2-1
        bits_down = res[i+1].log_cap - cpows[log_degree+1].log_cap
        res[idx+1] = mod_down_by(res[i+1], bits_down)
        res[idx+1] = mul(mk, res[idx+1], cpows[log_degree+1])
        res[idx+1] = rescale_by(res[idx+1], log_precision)
        idx += 1
    end
    res
end


function power_series_eager(
        mk::MultiplicationKey, cipher::Ciphertext, log_precision::Int,
        coeffs::Array{Float64, 1}, degree::Int)

    cpows = power_extended(mk, cipher, log_precision, degree)

    res = mul_by_const(cpows[1], coeffs[2], log_precision)
    res = add_const(res, coeffs[1])

    for i in 1:degree-1
        if abs(coeffs[i + 2]) > 1e-27 # TODO: why this limit?
            aixi = mul_by_const(cpows[i+1], coeffs[i + 2], log_precision)
            res = mod_down_to(res, aixi.log_cap)
            res = add(res, aixi)
        end
    end
    rescale_by(res, log_precision)
end


function power_series_lazy(
        mk::MultiplicationKey, cipher::Ciphertext, log_precision::Int,
        coeffs::Array{Float64, 1}, degree::Int)

    cpows = power_extended(mk, cipher, log_precision, degree)

    res = mul_by_const(cpows[1], coeffs[2], log_precision)
    res = add_const(res, coeffs[1])

    for i in 1:degree-1
        if abs(coeffs[i + 2]) > 1e-27 # TODO: why this limit?
            aixi = mul_by_const(cpows[i+1], coeffs[i + 2], log_precision)
            res = mod_down_by(res, res.log_cap - aixi.log_cap)
            res = add(res, aixi)
        end
    end

    res
end


power_series(mk, cipher, log_precision, coeffs, degree; lazy::Bool=false) =
    lazy ?
        power_series_lazy(mk, cipher, log_precision, coeffs, degree) :
        power_series_eager(mk, cipher, log_precision, coeffs, degree)


log_plus_one(mk, cipher, log_precision, degree; lazy::Bool=false) =
    power_series(
        mk, cipher, log_precision, [0, 1, -0.5, 1/3, -1/4, 1/5, -1/6, 1/7, -1/8, 1/9, -1/10],
        degree, lazy=lazy)


Base.exp(mk, cipher, log_precision, degree; lazy::Bool=false) =
    power_series(
        mk, cipher, log_precision,
        [1, 1, 0.5, 1/6, 1/24, 1/120, 1/720, 1/5040, 1/40320, 1/362880, 1/3628800],
        degree, lazy=lazy)


sigmoid(mk, cipher, log_precision, degree; lazy::Bool=false) =
    power_series(
        mk, cipher, log_precision,
        [1/2, 1/4, 0, -1/48, 0, 1/480, 0, -17/80640, 0, 31/1451520, 0],
        degree, lazy=lazy)
