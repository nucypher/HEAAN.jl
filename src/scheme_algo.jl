LOGARITHM = "Logarithm"
EXPONENT  = "Exponent" # exp(x)
SIGMOID = "Sigmoid" # sigmoid(x) = exp(x) / (1 + exp(x))


struct SchemeAlgo

    scheme :: Scheme
    taylorCoeffsMap :: Dict{String, Array{Float64, 1}}

    function SchemeAlgo(scheme::Scheme)

        tcmap = Dict{String, Array{Float64, 1}}()
        tcmap[LOGARITHM] = [0, 1, -0.5, 1/3, -1/4, 1/5, -1/6, 1/7, -1/8, 1/9, -1/10]
        tcmap[EXPONENT] = [1, 1, 0.5, 1/6, 1/24, 1/120, 1/720, 1/5040, 1/40320, 1/362880, 1/3628800]
        tcmap[SIGMOID] = [1/2, 1/4, 0, -1/48, 0, 1/480, 0, -17/80640, 0, 31/1451520, 0]

        new(scheme, tcmap)
    end
end


function powerOf2(algo::SchemeAlgo, cipher::Ciphertext, logp::Int, logDegree::Int)
    scheme = algo.scheme
    res = cipher

    for i in 0:logDegree-1
        res = square(scheme, res)
        res = reScaleBy(scheme, res, logp)
    end
    res
end


function powerOf2Extended(algo::SchemeAlgo, cipher::Ciphertext, logp::Int, logDegree::Int)
    scheme = algo.scheme
    res = Array{Ciphertext}(undef, logDegree + 1)
    res[0+1] = copy(cipher)
    for i in 1:logDegree
        res[i+1] = square(scheme, res[i-1+1])
        res[i+1] = reScaleBy(scheme, res[i+1], logp)
    end
    res
end


function power(algo::SchemeAlgo, cipher::Ciphertext, logp::Int, degree::Int)
    logDegree = floor(Int, log2(degree))
    po2Degree = 1 << logDegree

    scheme = algo.scheme

    res = powerOf2(algo, cipher, logp, logDegree)

    remDegree = degree - po2Degree
    if remDegree > 0
        tmp = power(algo, cipher, logp, remDegree)
        bitsDown = tmp.logq - res.logq
        tmp = modDownBy(scheme, tmp, bitsDown)
        res = mult(scheme, res, tmp)
        res = reScaleBy(scheme, res, logp)
    end
    res
end


function powerExtended(algo::SchemeAlgo, cipher::Ciphertext, logp::Int, degree::Int)
    scheme = algo.scheme
    logDegree = floor(Int, log2(degree))
    cpows = powerOf2Extended(algo, cipher, logp, logDegree)
    res = Array{Ciphertext}(undef, degree)
    idx = 0
    for i in 0:logDegree-1
        powi = (1 << i)
        res[idx+1] = copy(cpows[i+1])
        idx += 1
        for j in 0:powi-1-1
            bitsDown = res[j+1].logq - cpows[i+1].logq
            res[idx+1] = modDownBy(scheme, res[j+1], bitsDown)
            res[idx+1] = mult(scheme, res[idx+1], cpows[i+1])
            res[idx+1] = reScaleBy(scheme, res[idx+1], logp)
            idx += 1
        end
    end
    res[idx+1] = copy(cpows[logDegree+1])
    idx += 1
    degree2 = (1 << logDegree)
    for i in 0:degree-degree2-1
        bitsDown = res[i+1].logq - cpows[logDegree+1].logq
        res[idx+1] = modDownBy(scheme, res[i+1], bitsDown)
        res[idx+1] = mult(scheme, res[idx+1], cpows[logDegree+1])
        res[idx+1] = reScaleBy(scheme, res[idx+1], logp)
        idx += 1
    end
    res
end


function inverse(algo::SchemeAlgo, cipher::Ciphertext, logp::Int, steps::Int)
    scheme = algo.scheme
    cbar = negate(scheme, cipher)
    cbar = addConst(scheme, cbar, 1.0, logp)
    cpow = copy(cbar)
    tmp = addConst(scheme, cbar, 1.0, logp)
    tmp = modDownBy(scheme, tmp, logp)
    res = copy(tmp)
    for i in 1:steps-1
        cpow = square(scheme, cpow)
        cpow = reScaleBy(scheme, cpow, logp)
        tmp = copy(cpow)
        tmp = addConst(scheme, tmp, 1.0, logp)
        tmp = mult(scheme, tmp, res)
        tmp = reScaleBy(scheme, tmp, logp)
        res = copy(tmp)
    end
    res
end


function function_(algo::SchemeAlgo, cipher::Ciphertext, funcName, logp::Int, degree::Int)
    scheme = algo.scheme

    cpows = powerExtended(algo, cipher, logp, degree)

    dlogp = 2 * logp

    coeffs = algo.taylorCoeffsMap[funcName]

    res = multByConst(scheme, cpows[0+1], coeffs[1+1], logp)
    res = addConst(scheme, res, coeffs[0+1], dlogp)

    for i in 1:degree-1
        if abs(coeffs[i + 1 + 1]) > 1e-27 # TODO: why this limit?
            aixi = multByConst(scheme, cpows[i+1], coeffs[i + 1 + 1], logp)
            res = modDownTo(scheme, res, aixi.logq)
            res = add(scheme, res, aixi)
        end
    end
    reScaleBy(scheme, res, logp)
end


function functionLazy(algo::SchemeAlgo, cipher::Ciphertext, funcName, logp::Int, degree::Int)
    scheme = algo.scheme

    cpows = powerExtended(algo, cipher, logp, degree)

    dlogp = 2 * logp

    coeffs = algo.taylorCoeffsMap[funcName]

    res = multByConst(scheme, cpows[0+1], coeffs[1+1], logp)
    res = addConst(scheme, res, coeffs[0+1], dlogp)

    for i in 1:degree-1
        if abs(coeffs[i + 1 + 1]) > 1e-27 # TODO: why this limit?
            aixi = multByConst(scheme, cpows[i+1], coeffs[i + 1 + 1], logp)
            bitsDown = res.logq - aixi.logq
            res = modDownBy(scheme, res, bitsDown)
            res = add(scheme, res, aixi)
        end
    end

    res
end
