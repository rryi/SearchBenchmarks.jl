# script file to reproduce benchmarks in the analysis at
# TODO http ref
using SearchBenchmarks
using Test

path = "" # working directory including path delimiter at end
source = path*"7z.dll"

benchmark(200, 1000000, source, [naivesearch,bloom_v0], path*"(1) BIN julia impl much better than naive search",false)
benchmark(200, 1000000, source, [bloom_v0,bloom_v1,bloom_v2,bloom_v3], path*"(2) BIN julia bloom impl optimizations",false)
benchmark(200, 1100000, source, [bloom_v0,bloom_v3,bloom_v4], path*"(3) BIN do bloom test first",false)
benchmark(2000, 1000000, source, [naivesearch,bloom_v0,bloom_v3,bloom_v4], path*"(4) BIN degeneration on large patterns",false)
benchmark(2000, 1000000, source, [naivesearch,bloom_v0,bloom_v3,bloom_best], path*"(5) BIN limited bits in bloom filter",false)
benchmark(200, 1000000, source, [bloom_v0,bloom_v3,bloom_best,bmsearch], path*"(6) BIN boyer moore alternative",false)
benchmark(2000, 1000000, source, [naivesearch,bloom_v0,bloom_v3,bloom_best,bmsearch], path*"(7) BIN boyer moore on large patterns",false)
