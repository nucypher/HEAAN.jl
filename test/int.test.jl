using HEAAN: all_bits_mask, modulus, high_bit_mask, is_negative, num_bits, right_shift_rounded


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

    @testcase "right shift rounding" begin

        function ref_shift(x, shift, log_modulus)
            res = round(Int, x / 2^shift)
            # Since our integer range is (-q/2, q/2], not [-q/2, q/2],
            # there is a little kink if the result is -q/2.
            # We just add 1 in this case.
            if res == -(1 << (log_modulus - shift - 1))
                res += 1
            end
            res
        end

        log_modulus = 10
        shift = 5
        tp1 = BinModuloInt{Int, log_modulus}
        tp2 = BinModuloInt{Int, log_modulus - shift}

        for x in -2^(log_modulus-1)+1:2^(log_modulus-1)
            x_tp = convert(tp1, x)
            res = right_shift_rounded(x_tp, shift)
            ref = convert(tp2, ref_shift(x, shift, log_modulus))
            if res != ref
                @test_fail "Shift of $x failed: got $res, expected $ref"
                return
            end
        end
    end
end
