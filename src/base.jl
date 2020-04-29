# comparing _searchindex from base.search.jl with variants
# Author Robert Rudolph, Rudolph Consulting, Ahrensburg, Germany
import Base.ByteArray, Base._search_bloom_mask, Base._nthbyte


"""
data types supported for pattern and serarchable data.

We are investigating with respect to speed, so we restrict
to data structures with efficient byte access via position
and intentionally do not support AbstractVector and AbstractString
"""
const SearchSequence = Union{String,Vector{UInt8},Vector{Int8}}

#=
"structure used by many bloom-search-variants"
struct BloomSearch{T<:SearchSequence} <: SearchStructure
    pattern::T
    bloom_mask::UInt64
    bloom_skip::Int
    skip::Int # skip if last byte matches but another byte not
end
=#


@enum StatsFields begin
    SFinit=1        # time (ns) for initializing the search structure
    SFkernel        # time (ns) of search kernel execution
    SFloops         # no. of main loop runs of the search kernel
    SFtests         # no. of performed higher order (eg bloom) tests.
    SFskips         # no. of performed higher order (eg bloom) skips
    SFbits          # no. of bits set in the bloom filter
end


"Benchmark measurements. nothing ==> compiler strips measurement code"
const MaybeVector = Union{Nothing,Vector{Int}}


function bitcount(bitset::UInt64)
    count = 0
    while  bitset>0
        count += bitset&(1%UInt64)
        bitset = bitset >>>1
    end
    count
end


"benchmark measurement recording"
mutable struct Benchmark
    dict:: IdDict{Function,Array{Int,2}} # measurements
    size:: Int # no. of measurements per function in benchmark
    source:: String ## either UInt8,Int (alphabet end and seed) or a filename
    result::String
    function Benchmark(size::Int, source::String, result::String)
        d = IdDict{Function,Array{Int,2}}()
        new(d,size,source,result)
    end
end


# write benchmark measurements as CSV file
function Base.write(io::IO,s::Benchmark)
    for p in s.dict
        for j in Int(typemin(StatsFields)):Int(typemax(StatsFields))
            print(io,string(p.first),",",StatsFields(j))
            for v in p.second[:,j]
                print(io,", ",v)
            end
            println(io)
        end
    end
end

const nullmat = zeros(Int,0)

"record a line of measurements in the benchmark"
function record(stats::Benchmark, f::Function, i::Int, sv::MaybeVector)
    mat = get(stats.dict,f,nullmat)
    if mat===nullmat
        mat = zeros(Int,Int(stats.size),Int(typemax(StatsFields)))
        push!(stats.dict,f=>mat)
    end
    mat[i,:] = sv
end


"return stats series for function and field."
function line(stats::Benchmark,f::Function, field::StatsFields)
    m = get(stats.dict,f,nullmat)
    m[:,Int(field)]
end
