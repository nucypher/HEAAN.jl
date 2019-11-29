using HEAAN: all_bits_mask, modulus, high_bit_mask, is_negative, num_bits


@testgroup "Integer helper functions" begin

    @testcase "all_bits_mask()" begin
        @test all_bits_mask(Int, 10) == 2^10-1
    end

    @testcase "modulus()" begin
        @test modulus(Int, 10) == 2^10
    end

    @testcase "high_bit_mask()" begin
        @test high_bit_mask(Int, 10) == 2^9
    end

    @testcase "is_negative()" begin
        @test is_negative(2^15, 16) == false
        @test is_negative(2^15+1, 16) == true
        @test is_negative(2^16-1, 16) == true
    end

    @testcase "num_bits()" begin
        @test num_bits(UInt32(0)) == 0
        @test num_bits(UInt32(2^31)) == 32

        @test num_bits(Int32(0)) == 0
        @test num_bits(Int32(-2^31+1)) == 31
        @test num_bits(Int32(2^31-1)) == 31
        @test_throws DomainError num_bits(Int32(typemin(Int32)))

        @test num_bits(big(0)) == 0
        @test num_bits(big(1) << 32) == 33
        @test num_bits(big(1) << 234) == 235
        @test num_bits(-big(1) << 234) == 235
    end
end
