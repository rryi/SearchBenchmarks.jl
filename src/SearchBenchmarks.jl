module SearchBenchmarks

using Plots
import Random

include("base.jl") # basic types and utilities for benchmarking

# it follows the list of all search implementations.
# each file defines a search function with two methods:
# search initialization and search kernel (using initialization result).
include("bloom_v0.jl")
include("bloom_v1.jl")
include("bloom_v2.jl")
include("bloom_v3.jl")
include("bloom_v4.jl")
export bloom_v0, bloom_v1, bloom_v2, bloom_v3, bloom_v4
include("bloom_best.jl")
export bloom_best
# @inbounds versions to see effect of lowlevel optimizations
include("bloom_v0i.jl")
include("bloom_besti.jl")
export bloom_v0i, bloom_besti
include("bmsearch.jl")
export bmsearch
include("naivesearch.jl")
export naivesearch

# add more implementations...

include("benchmark.jl") # do and visualize benchmark

export benchmark, Benchmark
end # module
