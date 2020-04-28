"""
v1 but main loop until w-1 and additional end test
"""
function bloom_v2 end

function bloom_v2(t::SearchSequence)
    n = sizeof(t)
    bloom_mask = UInt64(0)
    skip = n - 1
    tlast = _nthbyte(t,n)
    for j in 1:n
        bloom_mask |= _search_bloom_mask(_nthbyte(t,j))
        if _nthbyte(t,j) == tlast && j < n
            skip = n - j -1
        end
    end
    return t,bloom_mask,bloom_skip,skip
end

function bloom_v2(s::SearchSequence, t::SearchSequence, i::Integer,sv::MaybeVector=nothing)
    bloom_v2(s,bloom_v2(t),s,i,sv)
end

function bloom_v2(s::SearchSequence,p::Tuple, i::Integer,sv::MaybeVector=nothing)
    DOSTATS = sv isa Nothing
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

    skip = n
    tlast = _nthbyte(t,n)
    bloom_mask = UInt64(_search_bloom_mask(tlast))
    for j in 1:n-1
        bloom_mask |= _search_bloom_mask(_nthbyte(t,j))
        if _nthbyte(t,j) == tlast
            skip = n - j
        end
    end
    if DOSTATS loops = 0 end
    if DOSTATS bloomtests = 0 end
    if DOSTATS bloomskips = 0 end
    i -= 1
    while i < w
        if DOSTATS loops += 1 end
        if _nthbyte(s,i+n) == tlast
            # check candidate
            j = 1
            while j < n
                if _nthbyte(s,i+j) != _nthbyte(t,j)
                    break
                end
                j += 1
            end
            # match found?
            if j == n
                if DOSTATS sv[Int(SFloops)] = loops; sv[Int(SFtests)] = bloomtests; sv[Int(SFskips)] = bloomskips; sv[Int(SFbits)] = bitcount(bloom_mask) end
                return i+1
            end
            # no match: skip and test bloom
            i += skip
            if DOSTATS bloomtests += 1 end
            if i<w && bloom_mask & _search_bloom_mask(_nthbyte(s,i+1+n)) == 0
                if DOSTATS bloomskips += 1 end
                i += n
            end
        else
            if DOSTATS bloomtests += 1 end
            # i<w is guaranteed here, comparison eliminated
            if bloom_mask & _search_bloom_mask(_nthbyte(s,i+1+n)) == 0
                if DOSTATS bloomskips += 1 end
                i += n
            end
        end
        i += 1
    end
    if i==w
        # test end match
        j = 1
        while j <= n
            if _nthbyte(s,i+j) != _nthbyte(t,j)
                break # not found
            end
            j += 1
        end
        if j > n
            if DOSTATS sv[Int(SFloops)] = loops; sv[Int(SFtests)] = bloomtests; sv[Int(SFskips)] = bloomskips; sv[Int(SFbits)] = bitcount(bloom_mask) end
            return i+1
        end
    end
    if DOSTATS sv[Int(SFloops)] = loops end
    if DOSTATS sv[Int(SFtests)] = bloomtests end
    if DOSTATS sv[Int(SFskips)] = bloomskips end
    if DOSTATS sv[Int(SFbits)] = bitcount(bloom_mask) end
    0
end
