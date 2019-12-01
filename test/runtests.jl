using Jute

include("int.test.jl")
include("float.test.jl")
include("rns.test.jl")
include("embedding.test.jl")

exit(runtests())
