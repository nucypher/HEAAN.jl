const logN = 16
const logQ = 600

const sigma = 3.2
const h = 64
const pbnd = 59
const kbar = 60
const kbar2 = 120
const logNh = logN - 1
const logQQ = 2 * logQ
const N = 1 << logN
const Nh = 1 << logNh
const M = N << 1
const nprimes = (2 + logN + 4 * logQ + pbnd - 1) ÷ pbnd
const Nnprimes = nprimes << logN
# const cbnd = (logQQ + NTL_ZZ_NBITS - 1) ÷ NTL_ZZ_NBITS # ???
const bignum = Int(0xfffffff)
const Q = BigInt(1) << logQ
const QQ = BigInt(1) << logQQ


function check_range(x::Array{BigInt, 1}, logq::Int)
    @assert !any(signbit.(x))
    nbs = maximum(num_bits.(x))
    @assert nbs <= logq
    #@assert nbs == logq "max nbs = $nbs, with logq = $logq"
end


function is_negative(x::BigInt, logq::Int)
    @assert !signbit(x)
    bit(x, logq - 1)
end

