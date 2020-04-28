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
    return t,bloom_mask,bloom_skip,skip
end

function bloom_v4(s::SearchSequence, t::SearchSequence, i::Integer,sv::MaybeVector=nothing)
    bloom_v4(s,bloom_v4(t),i,sv)
end

function bloom_v4(s::SearchSequence,p::Tuple,i::Integer,sv::MaybeVector=nothing)
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

    (t,bloom_mask,blook_skip,skip) = p

    if DOSTATS loops = 0 end
    if DOSTATS bloomtests = 0 end
    if DOSTATS bloomskips = 0 end
    i +=n-1
    while i <= m
        if DOSTATS loops += 1 end
        if DOSTATS bloomtests += 1 end
        if bloom_mask & _search_bloom_mask(_nthbyte(s,i)) == 0
            if DOSTATS bloomskips += 1 end
            i += n
        else
            if _nthbyte(s,i) == tlast
                # check candidate
                j = 1
                while _nthbyte(s,i-n+j) == _nthbyte(t,j)
                    j += 1
                    # match found?
                    if j == n
                        if DOSTATS sv[Int(SFloops)] = loops; sv[Int(SFtests)] = bloomtests; sv[Int(SFskips)] = bloomskips; sv[Int(SFbits)] = bitcount(bloom_mask) end
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
    if DOSTATS sv[Int(SFloops)] = loops end
    if DOSTATS sv[Int(SFtests)] = bloomtests end
    if DOSTATS sv[Int(SFskips)] = bloomskips end
    if DOSTATS sv[Int(SFbits)] = bitcount(bloom_mask) end
    0
end
