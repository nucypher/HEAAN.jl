struct Ciphertext

    ax :: Array{BigInt, 1}
    bx :: Array{BigInt, 1}

    logp :: Int
    logq :: Int

    n :: Int

    function Ciphertext(logp::Int, logq::Int, n::Int)
        ax = Array{BigInt}(undef, N)
        bx = Array{BigInt}(undef, N)
        new(ax, bx, logp, logq, n)
    end
end
