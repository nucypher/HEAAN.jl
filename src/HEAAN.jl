module HEAAN

using DarkIntegers
using Primes
using Random

include("int.jl")
include("float.jl")
include("compat.jl")
include("rns.jl")
include("embedding.jl")

include("params.jl")
include("polynomial.jl")
include("secret_key.jl")
include("public_key_set.jl")
include("ciphertext.jl")
include("operators.jl")
include("algo.jl")
include("bootstrap.jl")

export Params
export Context
export SecretKey
export encrypt
export decrypt
export EncryptionKey
export MultiplicationKey
export LeftRotationKey
export ConjugationKey
export add
export add_const
export mul
export imul
export mul_by_const
export mul_by_const_vec

export power_of_2
export power
export log_plus_one
export sigmoid

export BootstrapKey
export bootstrap

end
