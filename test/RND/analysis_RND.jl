# script file to reproduce benchmarks in the analysis at
# TODO http ref
using SearchBenchmarks
using Test

path = "" # working directory including path delimiter at end
source = "0x7f,123"

benchmark(200, 1000000, source, [naivesearch,bloom_v0], path*"(1) RND julia impl much better than naive search",false)
benchmark(200, 1000000, source, [bloom_v0,bloom_v1,bloom_v2,bloom_v3], path*"(2) RND julia bloom impl optimizations",false)
benchmark(200, 1000000, source, [bloom_v0,bloom_v3,bloom_v4], path*"(3) RND do bloom test first",false)
benchmark(2000, 1000000, source, [naivesearch,bloom_v0,bloom_v3,bloom_v4], path*"(4) RND degeneration on large patterns",false)
benchmark(2000, 1000000, source, [naivesearch,bloom_v0,bloom_v3,bloom_best], path*"(5) RND limited bits in bloom filter",false)
benchmark(200, 1000000, source, [bloom_v0,bloom_v3,bloom_best,bmsearch], path*"(6) RND boyer moore alternative",false)
benchmark(2000, 1000000, source, [naivesearch,bloom_v0,bloom_v3,bloom_best,bmsearch], path*"(7) RND boyer moore on large patterns",false)
