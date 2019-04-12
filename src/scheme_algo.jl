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
