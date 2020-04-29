# script file to reproduce benchmarks in the analysis at
# TODO http ref
using SearchBenchmarks
using Test

source = "somejuliasource.jl"
benchmark(100, 100000, source, [naivesearch,bloom_v0], "(1) TXT julia impl much better than naive search",false)
benchmark(100, 100000, source, [bloom_v0,bloom_v1,bloom_v2,bloom_v3], "(2) TXT julia bloom impl optimizations",false)
benchmark(100, 100000, source, [bloom_v0,bloom_v3,bloom_v4], "(3) TXT do bloom test first",false)
benchmark(1000, 100000, source, [bloom_v0,bloom_v3,bloom_v4], "(4) TXT degeneration on large patterns",false)
benchmark(1000, 100000, source, [bloom_v0,bloom_v3,bloom_best], "(5) TXT best-of-breed bloom",false)
benchmark(100, 100000, source, [bloom_v3,bloom_best,bmsearch], "(6) TXT boyer moore alternative",false)
benchmark(1000, 100000, source, [bloom_v3,bloom_best,bmsearch], "(7) TXT boyer moore on large patterns",false)
