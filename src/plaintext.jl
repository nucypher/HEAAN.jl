struct Plaintext

    mx :: Array{BigInt, 1}
    logp :: Int
    logq :: Int
    n :: Int

    function Plaintext(logp::Int, logq::Int, n::Int)
        mx = Array{BigInt}(undef, N)
        new(mx, logp, logq, n)
    end
end
