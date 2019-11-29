using Random
using HEAAN: float_to_integer, integer_to_float


function shift_bigfloat(x::BigFloat, shift::Int)
    # TODO: in NTL it's `MakeRR(x.x, x.e + shift)`,
    # but Julia does not give access to the significand.
    # Kind of hacky, but will work for now.
    xc = copy(x)
    xc.exp += shift
    xc
end


function float_to_integer_reference(x::Float64, shift::Int, log_modulus::Int)
    r = BigFloat(x)
    xp = shift_bigfloat(r, shift)
    res = round(BigInt, xp)
    signbit(res) ? (one(BigInt) << log_modulus) + res : res
end


function integer_to_float_reference(x::BigInt, shift::Int, log_modulus::Int)
    if x >= one(BigInt) << (log_modulus - 1)
        x -= one(BigInt) << log_modulus
    end
    xp = BigFloat(x)
    xp = shift_bigfloat(xp, -shift)
    convert(Float64, xp)
end


@testgroup "Float <-> integer conversions" begin


@testcase "Float to integer" begin
    rng = MersenneTwister(123)

    log_modulus = 150

    for i in 1:10000
        x = randn(rng) * 100

        shift = 20

        ref = float_to_integer_reference(x, shift, log_modulus)
        res = float_to_integer(BigInt, x, shift, log_modulus)

        if ref != res
            @test_fail "Converting $x, got $res, expected $ref"
            return
        end
    end
end


@testcase "Float to integer, corner cases" begin

    # zero
    x = float_to_integer(Int, 0.0, 15, 16)
    @test x == 0

    # biggest positive number
    x = float_to_integer(Int, 1.0, 15, 16)
    @test x == 2^15

    # smallest negative number
    x = float_to_integer(Int, -0.999969482421875, 15, 16)
    @test x == 2^15+1

    # biggest negative number
    x = float_to_integer(Int, -3.0517578125e-5, 15, 16)
    @test x == 2^16-1

    # very small number, rounded to zero
    x = float_to_integer(Int, 1e-10, 15, 16)
    @test x == 0

    # values outside the range
    @test_throws DomainError float_to_integer(Int, 1.01, 15, 16)
    @test_throws DomainError float_to_integer(Int, -1.0, 15, 16)

end


@testcase "Integer to float" begin
    rng = MersenneTwister(123)
    log_modulus = 100

    for i in 1:10000
        x = rand(rng, zero(BigInt):(one(BigInt)<<log_modulus)-one(BigInt))

        shift = 20

        ref = integer_to_float_reference(x, shift, log_modulus)
        res = integer_to_float(Float64, x, shift, log_modulus)

        if ref != res
            @test_fail "Converting $x, got $res, expected $ref"
            return
        end
    end
end


@testcase "Integer to float, corner cases" begin

    # zero
    x = integer_to_float(Float64, 0, 15, 16)
    @test x == 0.0

    # biggest positive number
    x = integer_to_float(Float64, 2^15, 15, 16)
    @test x == 1.0

    # smallest negative number
    x = integer_to_float(Float64, 2^15+1, 15, 16)
    @test x == -0.999969482421875

    # biggest negative number
    x = integer_to_float(Float64, 2^16-1, 15, 16)
    @test x == -3.0517578125e-5

    # very small number, rounded to zero
    x = integer_to_float(Float64, 1, 1040, 16)
    @test x == 0.0

    # values outside [0, modulus)
    @test_throws DomainError integer_to_float(Float64, 2^16, 15, 16)
    @test_throws DomainError integer_to_float(Float64, -1, 15, 16)

end


end
