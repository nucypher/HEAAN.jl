using Random
using HEAAN: EmbeddingPlan, embed, unembed


apply_polynomial(p, x) = sum(p .* x .^ collect(0:length(p)-1))


@testgroup "Embedding" begin


@testcase "forward" begin
    rng = MersenneTwister(123)
    log_max_len = 8

    plan = EmbeddingPlan(2^log_max_len)

    for log_len in 0:log_max_len
        len = 2^log_len
        a = randn(rng, len) + im * randn(rng, len)

        res = embed(plan, a)

        # `mod(, 4 * len)` is not necessary mathematically, but without it
        # `isapprox` fails because too much error is introduced at higher powers.
        ref = [apply_polynomial(a, exp(im * pi / (2 * len) * mod(5^j, 4 * len))) for j in 0:len-1]

        if !isapprox(res, ref)
            @test_fail "Embedding contract failed for length=$len"
            return
        end
    end
end


@testcase "inverse" begin
    rng = MersenneTwister(123)
    log_max_len = 8

    plan = EmbeddingPlan(2^log_max_len)

    for log_len in 0:log_max_len
        len = 2^log_len
        a = randn(rng, len) + im * randn(rng, len)

        res = unembed(plan, a)

        # `mod(, 4 * len)` is not necessary mathematically, but without it
        # `isapprox` fails because too much error is introduced at higher powers.
        ref = [apply_polynomial(res, exp(im * pi / (2 * len) * mod(5^j, 4 * len))) for j in 0:len-1]

        if !isapprox(a, ref)
            @test_fail "Unembedding contract failed for length=$len"
            return
        end
    end
end


end
