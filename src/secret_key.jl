struct SecretKey

    sx :: Array{BigInt, 1}

    function SecretKey(rng::MyRNG, ring::Ring)
        #sx = Array{BigInt}(undef, N)
        sx = sampleHWT(rng, ring)
        new(sx)
    end
end
