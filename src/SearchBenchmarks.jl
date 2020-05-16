module SearchBenchmarks

using Plots
import Random

include("base.jl") # basic types and utilities for benchmarking
export all_search_functions, Benchmark

# it follows the list of all search implementations.
# each file defines a search function with two methods:
# search initialization and search kernel (using initialization result).
include("naivesearch.jl")
export naivesearch
include("bloom_v0.jl")
include("bloom_v0i.jl")
include("bloom_v1.jl")
include("bloom_v2.jl")
include("bloom_v3.jl")
include("bloom_v4.jl")
export bloom_v0, bloom_v0i, bloom_v1, bloom_v2, bloom_v3, bloom_v4
include("bloom_best.jl")
include("bloom_best2.jl")
include("bloom_best3.jl")
include("bloom_bestc.jl")
export bloom_best, bloom_best2, bloom_best3,bloom_bestc
include("bmsearch.jl")
export bmsearch
#add more search implementations...

include("benchmark.jl") # do and visualize benchmark
export benchmark

end # module
