struct BootContext

    rpvec :: Array{Array{UInt64, 1}, 1}
    rpvecInv :: Array{Array{UInt64, 1}, 1}
    rp1 :: Array{UInt64, 1}
    rp2 :: Array{UInt64, 1}

    bndvec :: Array{Int, 1}
    bndvecInv :: Array{Int, 1}
    bnd1 :: Int
    bnd2 :: Int

    logp :: Int

    function BootContext(rpvec, rpvecInv, rp1, rp2, bndvec, bndvecInv, bnd1, bnd2, logp)
        new(rpvec, rpvecInv, rp1, rp2, bndvec, bndvecInv, bnd1, bnd2, logp)
    end
end
