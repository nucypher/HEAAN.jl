module HEAAN

using DarkIntegers
using Primes
using Random

include("int.jl")
include("float.jl")
include("rns.jl")
include("embedding.jl")
include("polynomial.jl")
include("random.jl")

include("params.jl")
export Params

include("secret_key.jl")
export SecretKey

include("public_keys.jl")
export EncryptionKey
export MultiplicationKey
export LeftRotationKey
export ConjugationKey

include("ciphertext.jl")
export encrypt
export decrypt

include("operators.jl")
export add
export sub
export add_const
export mul
export square
export negate
export imul
export mul_by_const
export mul_by_const_vec
export div_by_po2
export mod_down_by
export rescale_by

include("algo.jl")
export power
export sigmoid
export log_plus_one
export power_series

include("bootstrap.jl")
export BootstrapKey
export bootstrap

end
