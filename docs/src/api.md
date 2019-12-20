# API reference

```@meta
CurrentModule = HEAAN
```


## Scheme parameters

```@docs
Params
```


## Keys

```@docs
SecretKey
EncryptionKey
MultiplicationKey
LeftRotationKey
ConjugationKey
```


## Encryption/decryption

```@docs
Ciphertext
encrypt
decrypt
```


## Arithmetic operations

```@docs
add
add_const
sub
mul
square
negate
imul
mul_by_const
mul_by_const_vec
div_by_po2
circshift
conj
inv
power
power_series
sigmoid
log_plus_one
exp
```


## Computation resource management

```@docs
mod_down_by
mod_down_to
rescale_by
BootstrapKey
bootstrap
```
