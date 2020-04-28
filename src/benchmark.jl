

"""
run benchmark measurements
"""
function runbench(f::Function,s,t, b::Benchmark, stats::Bool)
    #print(string(f),": ")
    sf = zeros(Int,Int(typemax(StatsFields)))
    starttime = time_ns()
    pattern =f(t)
    inittime =  time_ns()
    if stats
        pos = f(s,pattern,1,sf)
    else
        pos = f(s,pattern,1,nothing)
    end
    kerneltime = time_ns()
    sf[Int(SFinit)] = (inittime - starttime)%Int
    sf[Int(SFkernel)] = (kerneltime-inittime)%Int
    record(stats,f,sizeof(t),sf)
    pos
end

#=
"""
    benchmark(patternsize::Int, textsize::Int, source::String, functions::Vector{Function}, result::String, stats::Bool)

run a benchmark: call search functions for specified pattern and sequence to search,
and plot timing results. Plots are saved as graphic files, filename is made up
from the parameters

patternsize: no. of bytes for the longest pattern.
patterns will vary from size 3 up to patternsize

textsize: number of bytes in the text to search in.

source: either a filename or an UInt8 and an Int separated by comma.
If source begins with 0x, it must consist of an UInt8 and an Int, separated
by a comma, in julia literal encoding. e.g.: "0x7f,1".
benchmark will generate a random sequence of bytes of length textsize, and
the UInt8 specifies the highest byte value to use.
The Int specifies a seed for the random generator, to make runs reproducible.

functions: the names of all search functions to benchmark.
Proper function definitions must be included in the code.

result: filename with complete path but no file extension. 3 files
are written: "$result.csv" contains timing and satistics in a CSV file,
"$result.png" the kernel timing plot, "$result.init.png" the initialitaion
timings plot.

stats: if true, statistics will be collected
and saved in a file with this name, in CSV format.
"""
=#
function benchmarkx(patternsize::Int, textsize::Int, source::String, functions::Vector{Function}, result::String, stats::Bool)
    stats = Benchmark(textsize,source,result)
    s = Vector{UInt8}(undef,textsize)
    if startswith(source,"0x")
        # generate random sequence
        args = split(source,',')
        alphabetsize = parse(UInt8,args[1])
        seed = parse(Int,args[2])
        Random.seed!(seed)
        rand!(s,0x00:alphabetsize) # array to search in
    else
        # read file
        f = open(source)
        read!(f,s)
        close(f)
    end
    0
end


function benchmark(patternsize::Int, textsize::Int, source::String, functions::Vector{Function}, result::String, stats::Bool)
    stats = Benchmark(textsize,source,result)
    s = Vector{UInt8}(undef,textsize)
    if startswith(source,"0x")
        # generate random sequence
        args = split(source,',')
        alphabetsize = parse(UInt8,args[1])
        seed = parse(Int,args[2])
        Random.seed!(seed)
        rand!(s,0x00:alphabetsize) # array to search in
    else
        # read file
        f = open(source)
        read!(f,s)
        close(f)
    end
    #ss = String(copy(s))
    for j in 4:patternsize
        tj = s[textsize-patternsize+1:textsize-patternsize+j]
        #println()
        #println("benchmark with pattern size $j")
        last = -1
        for f in functions
            pos = runbench(f,s,tj,b,stats)
            if pos!=last && last>=0
                @error "pos mismatch: f=$f, pos=$pos, last=$last"
            end
        end
    end
    if stats
        io = open(result*".csv", "w")
        print(io,b)
        close(io)
    end

    # plot
    ##pp =plot(line(stats,_searchindex_naive,SFtime), label=string(_searchindex_naive))
    plotinit =plot()
    title!("$result: initialization timing")
    #for f in [_searchindex_naive,_searchindex_julia,_searchindex_v1,_searchindex_v2,_searchindex_v3,_searchindex_v4,_searchindex_v5 ]
    for f in functions
        plot!(line(stats,f,SFinit), label=string(f))

    end
    savefig(plotinit,result*"_init")
    plotkernel =plot()
    title!("$result: kernel timing")
    #for f in [_searchindex_naive,_searchindex_julia,_searchindex_v1,_searchindex_v2,_searchindex_v3,_searchindex_v4,_searchindex_v5 ]
    for f in functions
        plot!(line(stats,f,SFkernel), label=string(f))

    end
    savefig(plot,result)
    display(plotinit)
    display(plotkernel)

end

#benchmark(128,3,100,"r:\\testfile.txt")
