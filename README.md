# A Julia implementation of the HEAAN scheme

Master branch: [![CircleCI](https://circleci.com/gh/nucypher/HEAAN.jl.svg?style=svg)](https://circleci.com/gh/nucypher/HEAAN.jl) [![codecov](https://codecov.io/gh/nucypher/HEAAN.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/nucypher/HEAAN.jl)


## Sources

This is an implementation of:

J. H. Cheon, A. Kim, M. Kim, and Y. Song, *"Homomorphic Encryption for Arithmetic of Approximate Numbers,"* [Lecture Notes in Computer Science, Vol. 10624 LNCS, pp. 409–437 (2017)](https://doi.org/10.1007/978-3-319-70694-8_15),

and partially

J. H. Cheon, K. Han, A. Kim, M. Kim, and Y. Song, *"Improved Bootstrapping for Approximate Homomorphic Encryption,"* [eprint 2018/1043](https://eprint.iacr.org/2018/1043).

For the further development of this scheme one can refer to

J. H. Cheon, K. Han, A. Kim, M. Kim, and Y. Song, *"A Full RNS Variant of Approximate Homomorphic Encryption,"* [eprint 2018/931](https://eprint.iacr.org/2018/931.pdf).

The implementation is based on the [reference C++ code](https://github.com/snucrypto/HEAAN).


## A simple example

The following example illustrates basics of working with the library.
For more examples see `examples/basic.jl` and `test/api.test.jl` in the repository.

First, we import the required libraries and initialize constants:

```julia
using Random
using HEAAN

# Polynomial length and full modulus determine the security of the scheme
params = Params(log_polynomial_length=8, log_lo_modulus=300)

# Vector size
# The maximum vector size is polynomial_length / 2
n = 2^6

# Initial precision for ciphertexts (that is, the absolute precision is 1/2^30)
log_precision = 30

# Initial precision cap for ciphertexts
# Gets consumed during e.g. multiplication
log_cap = 200

rng = MersenneTwister(123)
```

Then we need to create keys. A secret key is required to decrypt ciphertexts and initialize public keys.
There are several different public keys in HEAAN; in this example we will only use two: the one required for encryption, and the one required for multiplication.
Addition of ciphertexts can be performed without a key.

```julia
secret_key = SecretKey(rng, params)

# Public key used for encryption
enc_key = EncryptionKey(rng, secret_key)

# Public key used for multiplication
mul_key = MultiplicationKey(rng, secret_key)
```

We create two complex-valued vectors of length `n` and initialize the reference array.

```julia
v1 = rand(rng, n) + im * rand(rng, n)
v2 = rand(rng, n) + im * rand(rng, n)

# Reference calculation
ref = x .* y .+ x
```

Use `enc_key` to create initial ciphertexts.
`log_precision` specifies the absolute precision for encrypted values; `log_cap` is essentially the "resource" that gets spent during arithmetic operations on ciphertexts.
Naturally, `log_cap` should be greater than `log_precision`.

```julia
# Encrypt the initial vectors
c1 = encrypt(rng, enc_key, v1, log_precision, log_cap)
c2 = encrypt(rng, enc_key, v2, log_precision, log_cap)

println("Before: precision=$(c1.log_precision), cap=$(c1.log_cap)")

# output
Before: precision=30, cap=200
```

We start from performing the (elementwise) multiplication.

```julia
t1 = mul(mul_key, c1, c2)

println("After multiplication: precision=$(t1.log_precision), cap=$(t1.log_cap)")

# output
After multiplication: precision=60, cap=200
```

Note that the precision of the result is the sum of the precisions of ciphertexts, while the cap remained the same.
This means that you cannot just multiply indefinitely --- after several iterations you will lose the encrypted information.
This problem is solved by `bootstrap()` which increases the difference between `precision` and `log_cap`, essentially adding the "computation resource".

Another consequence of this is that now `t1` and `c1` which we want to add together have different `log_cap`.
In order to `add()` them, both their `log_cap` and `log_precision` must be equal.
There are two functions that help with that: `rescale_by()` and `mod_down_by()`.
The former decreases both `log_cap` and `log_precision` by the same amount; the latter just decreases `log_cap`.
As it happens, in our case we need both of them.

```julia
t2 = rescale_by(t1, 30)

println("After rescale: precision=$(t2.log_precision), cap=$(t2.log_cap)")

t3 = mod_down_by(c1, 30)

println("After mod_down: precision=$(t3.log_precision), cap=$(t3.log_cap)")

# output
After rescale: precision=30, cap=170
After mod_down: precision=30, cap=170
```

Now that the parameters are equalized, we can call `add()` and decrypt the result.

```julia
cres = add(t2, t3)

res = decrypt(secret_key, cres)
ref = reference(v1, v2)

for i in 1:n
    println("$i-th element: diff=$(abs(res[i] - ref[i]))")
end

# output
1-th element: diff=1.8071040828704262e-8
2-th element: diff=4.8677934407457675e-8
3-th element: diff=3.612811639903806e-8
4-th element: diff=2.9650269405023588e-8
...
```

As you can see, the calculation is accurate to about 25 bits.
