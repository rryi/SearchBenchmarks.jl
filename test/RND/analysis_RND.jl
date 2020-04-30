# script file to reproduce benchmarks in the analysis at
# TODO http ref
using SearchBenchmarks
using Test

source = "0x7f,17"
benchmark(100, 1000000, source, [naivesearch,bloom_v0], "(1) RND julia impl much better than naive search",false)
benchmark(100, 1000000, source, [bloom_v0,bloom_v1,bloom_v2,bloom_v3], "(2) RND julia bloom impl optimizations",false)
benchmark(100, 1000000, source, [bloom_v0,bloom_v3,bloom_v4], "(3) RND do bloom test first",false)
benchmark(1000, 1000000, source, [bloom_v0,bloom_v3,bloom_v4], "(4) RND degeneration on large patterns",false)
benchmark(1000, 1000000, source, [bloom_v0,bloom_v3,bloom_best], "(5) RND best-of-breed bloom",false)
benchmark(100, 1000000, source, [bloom_v3,bloom_best,bmsearch], "(6) RND boyer moore alternative",false)
benchmark(1000, 1000000, source, [bloom_v3,bloom_best,bmsearch], "(7) RND boyer moore on large patterns",false)
