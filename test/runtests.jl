using SearchBenchmarks
using Test
import Random

function module_functions(modname)
    list = Symbol[]
    for nm in names(modname)
        typeof(eval(nm)) == Function && push!(list,nm)
    end
    return list
end

@testset "SearchBenchmarks.jl" begin
    # test constants
    test_textsize = 10000 # length of searchable array
    test_highestbyte = 0x7f # highest value in alphabet (must be UInt8)
    test_patterns = 8 # no.of pattern sizes to test
    test_seed = 2345 # random generator seed
    # get search functions
    searches = all_search_functions()
    #@info "# search functions: ",searches

    # compare search results of all functions
    f0 = pop!(searches)
    s = Vector{UInt8}(undef,test_textsize)
    Random.seed!(test_seed)
    Random.rand!(s,0x00:test_highestbyte) # array to search in
    psize = 1
    for pi in 1:test_patterns
        #@info "pattern size", psize
        t = s[test_textsize-psize+1:test_textsize]
        pos0 = 0
        while true
            pos = pos0 + 1
            pos0 = f0(s,t,pos)
            for f in searches
                posf = f(s,t,pos)
                if pos0!=posf
                    # error! report params
                    @error "test fails on $f0($s,$t,$pos) == $f($s,$t,$pos): $pos0 == $posf"
                end
                @test pos0 == posf
                #@info "$f(s,t,pos) = $pp"
            end
            if pos0 == 0
                break
            end
        end
        psize = psize*3-1 # sequence 1 2 5 14 41 122 ...
    end
    # test benchmarking
    benchmark(5, 100, "0x7f,1234", searches, "SearchBenchmarks_test_result",false)
    benchmark(5, 100, "SearchBenchmarks_test_result.png", searches, "SearchBenchmarks_test_result2",true)
end
