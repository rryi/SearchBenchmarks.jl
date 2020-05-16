"""
variant: test bloomfilter first
"""
function bloom_v4 end


function bloom_v4(t::SearchSequence)
    n = sizeof(t)
    skip = n
    tlast = _nthbyte(t,n)
    bloom_mask = UInt64(_search_bloom_mask(tlast))
    for j in 1:n-1
        bloom_mask |= _search_bloom_mask(_nthbyte(t,j))
        if _nthbyte(t,j) == tlast
            skip = n - j
        end
    end
    return t,bloom_mask,skip,tlast
end

function bloom_v4(s::SearchSequence, t::SearchSequence, i::Integer)
    bloom_v4(s,bloom_v4(t),i)
end

function bloom_v4(s::SearchSequence,p::Tuple,i::Integer,sv::MaybeVector=nothing)
    (t,bloom_mask,skip,tlast) = p # replaces original code, commented out
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

    if DOSTATS loops = 0 end
    if DOSTATS bloomtests = 0 end
    if DOSTATS bloomskips = 0 end
    i +=n-1
    while i <= m
        if DOSTATS loops += 1 end
        if DOSTATS bloomtests += 1 end
        silast = _nthbyte(s,i)
        if bloom_mask & _search_bloom_mask(silast) == 0
            if DOSTATS bloomskips += 1 end
            i += n
        else
            if silast == tlast
                # check candidate
                j = 1
                while _nthbyte(s,i-n+j) == _nthbyte(t,j)
                    j += 1
                    # match found?
                    if j == n
                        if DOSTATS recordcase(sv, loops, bloomtests, bloomskips, bitcount(bloom_mask), skip) end
                        return i-n+1
                    end
                end
                # no match: skip and test bloom
                i += skip
            else
                i +=1
            end
        end
    end
    if DOSTATS recordcase(sv, loops, bloomtests, bloomskips, bitcount(bloom_mask), skip) end
    0
end
