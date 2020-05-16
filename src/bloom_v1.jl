"""
v0 with improved skip after a last byte match.

v0 does bloom test at current position +1 and uses skip only if
bloom test allows not for skipping.

v1 uses skip and tests bloom at current position + skip
"""


function bloom_v1(t::SearchSequence)
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

function bloom_v1(s::SearchSequence, t::SearchSequence, i::Integer)
    bloom_v1(s,bloom_v1(t),i)
end

function bloom_v1(s::SearchSequence,p::Tuple, i::Integer,sv::MaybeVector=nothing)
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
            # match found?
            if j == n - 1
                if DOSTATS recordcase(sv, loops, bloomtests, bloomskips, bitcount(bloom_mask), skip) end
                return i+1
            end

            # no match. skip forward and try bloom
            i += skip
        end
        if i < w
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
