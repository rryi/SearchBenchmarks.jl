

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
    if sizeof(t)>2
        sf[Int(SFinit)] = (inittime - starttime)%Int
        sf[Int(SFkernel)] = (kerneltime-inittime)%Int
        # patternsize 2 is 1st run, ignore it (may include compilation)
        record(b,f,sizeof(t),sf)
    end
end


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
        last = -1
        for f in functions
            pos = runbench(f,s,tj,b,stats)
            if pos!=last && last>=0
                @error "pos mismatch: f=$f, pos=$pos, last=$last"
            end
        end
    end
    io = open(result*".csv", "w")
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
    title!("$result: kernel timing")
    for f in functions
        plot!(line(b,f,SFkernel), label=string(f))
    end
    png(result)
    display(plotkernel)
end

#benchmark(128,3,100,"r:\\testfile.txt")
