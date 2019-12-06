using Jute

include("int.test.jl")
include("float.test.jl")
include("rns.test.jl")
include("embedding.test.jl")
include("api.test.jl")

exit(runtests())
