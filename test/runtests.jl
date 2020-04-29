using SearchBenchmarks
using Test

@testset "SearchBenchmarks.jl" begin
    # TODO

    benchmark(5, 1000, "0x7f,2000", [naivesearch,bloom_v0], "SearchBenchmarks_test_result",true)
    benchmark(5, 1000, "SearchBenchmarks_test_result.png", [naivesearch,bloom_v0,bloom_v1,bloom_v2,bloom_v3,bloom_v4,bloom_best,bmsearch], "SearchBenchmarks_test_result2",false)

end
