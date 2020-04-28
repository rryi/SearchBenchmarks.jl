module SearchBenchmarks

using Plots

include("base.jl") # basic types and utilities for benchmarking

# it follows the list of all search implementations.
# each file defines a search function with two methods:
# search initialization and search kernel (using initialization result).
include("bloom_v0.jl")
include("bloom_v1.jl")
include("bloom_v2.jl")
include("bloom_v3.jl")
include("bloom_v4.jl")
include("bloom_best.jl")
include("bmsearch.jl")
include("naivesearch.jl")

# add more implementations...

include("benchmark.jl") # do and visualize benchmark

end # module
