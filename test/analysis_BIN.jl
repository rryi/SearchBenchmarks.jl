# script file to reproduce benchmarks in the analysis at
# TODO http ref
using SearchBenchmarks
using Test

source = "7z.dll"
benchmark(100, 1000000, source, [naivesearch,bloom_v0], "(1) BIN julia impl much better than naive search",false)
benchmark(100, 1000000, source, [bloom_v0,bloom_v1,bloom_v2,bloom_v3], "(2) BIN julia bloom impl optimizations",false)
benchmark(100, 1000000, source, [bloom_v0,bloom_v3,bloom_v4], "(3) BIN do bloom test first",false)
benchmark(1000, 1000000, source, [bloom_v0,bloom_v3,bloom_v4], "(4) BIN degeneration on large patterns",false)
benchmark(1000, 1000000, source, [bloom_v0,bloom_v3,bloom_best], "(5) BIN best-of-breed bloom",false)
benchmark(100, 1000000, source, [bloom_v3,bloom_best,bmsearch], "(6) BIN boyer moore alternative",false)
benchmark(1000, 1000000, source, [bloom_v3,bloom_best,bmsearch], "(7) BIN boyer moore on large patterns",false)
