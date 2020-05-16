#=
general contract for search functions f in SearchBenchmarks:

    f(pattern::SearchSequence)::Tuple

initializes a search structure as Tuple, which includes pattern and
derived data for higher order search functions.
It is put in a separate method to allow for simpler timing
measurement of initialization code.

It is used in the search "kernel", which consists of the main
search loop:

    f(ss::Tuple,text::SearchSequence,start::Int)

=#


"""
This is the original implementation in julia Base,
method _searchindex(s::SearchSequence, t::SearchSequence, i::Integer).

Some technical changes are:
 * split into initialization and kernel
 * benchmark measurements, compiled away if last parameter is nothing
 * code is valid for String and ByteArray sequences
"""
function bloom_v0 end

function bloom_v0(t::SearchSequence)
    n = sizeof(t)
    bloom_mask = UInt64(0)
    skip = n - 1
    tlast = _nthbyte(t,n)
    for j in 1:n
        bloom_mask |= _search_bloom_mask(_nthbyte(t,j))
        if _nthbyte(t,j) == tlast && j < n
            skip = n - j - 1
        end
    end
    return t,bloom_mask,skip,tlast
end

function bloom_v0(s::SearchSequence, t::SearchSequence, i::Integer)
    bloom_v0(s,bloom_v0(t),i)
end

function bloom_v0(s::SearchSequence, p::Tuple,i::Integer,sv::MaybeVector=nothing)
    ( t,bloom_mask,skip,tlast) = p
    DOSTATS = !(sv isa Nothing)
    n = sizeof(t)
    m = sizeof(s)

    if n == 0
        return 1 <= i <= m+1 ? max(1, i) : 0
    elseif m == 0
        return 0
    elseif n == 1
        return something(findnext(isequal(_nthbyte(t,1)), s, i), 0)
    end

    w = m - n
    if w < 0 || i - 1 > w
        return 0
    end

#=
    bloom_mask = UInt64(0)
    skip = n - 1
    tlast = _nthbyte(t,n)
    for j in 1:n
        bloom_mask |= _search_bloom_mask(_nthbyte(t,j))
        if _nthbyte(t,j) == tlast && j < n
            skip = n - j - 1
        end
    end
=#
    if DOSTATS loops = 0 end
    if DOSTATS bloomtests = 0 end
    if DOSTATS bloomskips = 0 end
    i -= 1
    while i <= w
        if DOSTATS loops += 1 end
        if _nthbyte(s,i+n) == tlast
            # check candidate
            j = 0
            while j < n - 1
                if _nthbyte(s,i+j+1) != _nthbyte(t,j+1)
                    break
                end
                j += 1
            end

            # match found
            if j == n - 1
                if DOSTATS recordcase(sv, loops, bloomtests, bloomskips, bitcount(bloom_mask), skip) end
                return i+1
            end

            # no match, try to rule out the next character
            if DOSTATS bloomtests += 1 end
            if i < w && bloom_mask & _search_bloom_mask(_nthbyte(s,i+n+1)) == 0
                if DOSTATS bloomskips += 1 end
                i += n
            else
                i += skip
            end
        elseif i < w
            if DOSTATS bloomtests += 1 end
            if bloom_mask & _search_bloom_mask(_nthbyte(s,i+n+1)) == 0
                if DOSTATS bloomskips += 1 end
                i += n
            end
        end
        i += 1
    end
    if DOSTATS recordcase(sv, loops, bloomtests, bloomskips, bitcount(bloom_mask), skip) end
    0
end
