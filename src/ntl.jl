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


function bit(x::BigInt, i::Int)
    x & (one(BigInt) << i) != 0
end
