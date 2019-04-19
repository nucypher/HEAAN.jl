# Analogues of NTL functions

#=
function RandomBits_ZZ(len::Int)
    b = bitrand(len)
    reduce(|, (one(BigInt) .<< collect(0:len-1)) .* b)
end


function RandomBits_long(len::Int)
    b = bitrand(len)
    reduce(|, (one(Int) .<< collect(0:len-1)) .* b)
end
=#

# TODO: the NTL AddMod() should only be applied to numbers in range 0...q-1.
# In HEAAN, it is often applied to negative numbers, or numbers > q and somehow works,
# Even thougn in this case AddMod does not do anything, or just subtracts q, respectively.
# This function tries to imitate this behavior.
function AddMod(x::BigInt, y::BigInt, q::BigInt)
    res = x + y
    if res >= q
        #rem(res, q)
        res - q
    else
        res
    end
end


function MulMod(x::BigInt, y::BigInt, q::BigInt)
    mod(x * y, q)
end


function to_RR(x::Float64)
    BigFloat(x)
end


function to_RR(x::BigInt)
    BigFloat(x)
end


function MakeRR_shift(x::BigFloat, shift::Int)
    # TODO: in NTL it's `MakeRR(x.x, x.e + shift)`,
    # but Julia does not give access to the significand.
    # Kind of hacky, but will work for now.
    xc = copy(x)
    xc.exp += shift
    xc
end


function RoundToZZ(x::BigFloat)
    round(BigInt, x)
end


function bit(x::BigInt, i::Int)
    x & (one(BigInt) << i) != 0
end


function to_double(x::BigFloat)
    convert(Float64, x)
end


function NumBits(x::BigInt)
    # TODO: slow, but works
    if iszero(x)
        0
    elseif x < 0
        NumBits(-x)
    else
        ceil(Int, log2(x + one(BigInt)))
    end
end
