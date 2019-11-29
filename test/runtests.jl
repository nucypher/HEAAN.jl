using Jute

include("int.test.jl")
include("float.test.jl")
include("rns.test.jl")

exit(runtests())
