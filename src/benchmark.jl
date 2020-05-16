

"""
run benchmark measurements
"""
function runbench(f::Function,s,t, b::Benchmark, stats::Bool)
    #print(string(f),": ")
    sf = zeros(Int,Int(typemax(StatsFields)))
    starttime =  time_ns()
    pattern = f(t)
    inittime =  time_ns()
    pos = 0
    if stats
        while (pos = f(s,pattern,pos+1,sf))>0 end
    else
        while (pos = f(s,pattern,pos+1,nothing))>0 end
    end
    kerneltime = time_ns()
    if sizeof(t)>3 # do not report pattern size up to 3
        sf[Int(SFinit)] = (inittime - starttime)%Int
        sf[Int(SFkernel)] = (kerneltime-inittime)%Int
        # patternsize 2 is 1st run, ignore it (may include compilation)
        record(b,f,sizeof(t),sf)
    end
    nothing
end


"""
    benchmark(patternsize::Int, textsize::Int, source::String, functions::Vector{Function}, result::String, stats::Bool)

patternsize: largest size of search pattern in benchmark loop.
Loop runs from 2 (not measured) up to patternsize

textsize: length of the byte sequence to search in. Source must supply al least
so many bytes

source: either an UInt8 and an Int literal, separated by comma, or a file path.
If source starts with "0x" benchmark expects an UInt8 and Int litaral. The UInt8
specifies the highest value in the alphabet for a randomly generated search
sequence, the Int gives a random seed.

functions: the list of search functions to include in the benchmark

result: a file path. benchmarks writes files <result>.png and <result>.csv.
The former file shows the search kernel timings, the latter saves all collected
statistics in CSV format, including search timings.

stats: if true, additional measures are talen, like the number of main loop
cycles or the number of successful bloom skips. Thanks to multiple dispatch,
code collecting statistics is optimized away if stats==false. Because we are
benchmarking very lowlevel ans simple operations, code collecting statistics
can add overhead which results in senseless benchmark timings. Use stats==true
mostly in the design and anaysis phase of a search function, and switch stats
off for final benchmarking.

"""
Base.@propagate_inbounds function benchmark(patternsize::Int, textsize::Int, source::String, functions::Vector{Function}, result::String, stats::Bool)
    b = Benchmark(patternsize,source,result)
    s = Vector{UInt8}(undef,textsize)
    if startswith(source,"0x")
        # generate random sequence
        args = split(source,',')
        alphabetsize = parse(UInt8,args[1])
        seed = parse(Int,args[2])
        Random.seed!(seed)
        Random.rand!(s,0x00:alphabetsize) # array to search in
    else
        # read file
        f = open(source)
        read!(f,s)
        close(f)
    end
    #ss = String(copy(s))
    for j in 2:patternsize
        tj = s[textsize-patternsize+1:textsize-patternsize+j]
        #println()
        #println("benchmark with pattern size $j")
        for f in functions
            runbench(f,s,tj,b,stats)
        end
    end
    io = open(result * ".csv", "w")
    write(io,b)
    close(io)

    #= initializing timings need more repititions to be statistically valid
    plotinit =plot()
    title!("$result: initialization timing")
    for f in functions
        plot!(line(b,f,SFinit), label=string(f))
    end
    png(result*"_init")
    =#
    plotkernel =plot()
    # remove (windows) path part from result for title to keep plot title intact
    title = result
    endpath = findlast('\\',result)
    if endpath != nothing
        title = result[endpath+1:end]
    end
    title!(title*": kernel timing")
    for f in functions
        plot!(line(b,f,SFkernel), label=string(f))
    end
    png(result)
    display(plotkernel)
    nothing
end

#benchmark(128,3,100,"r:\\testfile.txt")
