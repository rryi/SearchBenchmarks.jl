# script file to reproduce benchmarks in the analysis at
# TODO http ref
using SearchBenchmarks
using Test

path = "" # working directory including path delimiter at end
source = path*"somejuliasource.jl"

benchmark(200, 100000, source, [naivesearch,bloom_v0], path*"(1) TXT julia impl much better than naive search",false)
benchmark(200, 100000, source, [bloom_v0,bloom_v1,bloom_v2,bloom_v3], path*"(2) TXT julia bloom impl optimizations",false)
benchmark(200, 100000, source, [bloom_v0,bloom_v3,bloom_v4], path*"(3) TXT do bloom test first",false)
benchmark(2000, 100000, source, [naivesearch,bloom_v0,bloom_v3,bloom_v4], path*"(4) TXT degeneration on large patterns",false)
benchmark(2000, 100000, source, [naivesearch,bloom_v3,bloom_v4,bloom_best], path*"(5) TXT limited bits in bloom filter",false)
benchmark(200, 100000, source, [bloom_v0,bloom_v3,bloom_best,bmsearch], path*"(6) TXT boyer moore alternative",false)
benchmark(2000, 100000, source, [naivesearch,bloom_v0,bloom_v3,bloom_best,bmsearch], path*"(7) TXT boyer moore on large patterns",false)
