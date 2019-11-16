randomReal(rng::MyRNG, bound::Float64 = 1.0) = myrand_float(rng) * bound
randomComplex(rng::MyRNG, bound::Float64 = 1.0) =
    randomReal(rng, bound) + im * randomReal(rng, bound)
randomComplexArray(rng::MyRNG, n::Int, bound::Float64 = 1.0) =
    [randomComplex(rng, bound) for i in 1:n]

function randomCircle(rng::MyRNG, anglebound::Float64 = 1.0)
    angle_ = randomReal(rng, anglebound)
    M_PI = Float64(pi)
    mycos(angle_ * 2 * M_PI) + im * mysin(angle_ * 2 * M_PI)
end


float_to_integer(x::Union{BigFloat, Float64}, shift::Int) =
    float_to_integer(BigInt, x, shift)
