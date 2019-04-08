struct Key

    rax :: Array{UInt64, 1}
    rbx :: Array{UInt64, 1}

    function Key()
        rax = Array{UInt64}(undef, Nnprimes)
        rbx = Array{UInt64}(undef, Nnprimes)
        new(rax, rbx)
    end
end
