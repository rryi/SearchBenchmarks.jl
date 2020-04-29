"""
Optimized search function, using benchmark insights:

    * use "bloom first" (bloom_v4) for small patterns

    * restrict skips on large patterns, avoiding filter degeneration
      and also reducing initialization effort

"""
function bloom_besti end

function bloom_besti(t::SearchSequence)
    n = sizeof(t)
    skip = n
    tlast = _nthbyte(t,n)
    bloom_mask = UInt64(_search_bloom_mask(tlast))
    bloom_bits = 1 # no. of bits set in bloom filter if <= 32
    bloom_skip = 1 # no. of bytes to skip if byte not in filter
    j = n
    while (j-=1)>=1  # j=n is already processed
        tj = _nthbyte(t,j)
        if tj == tlast && skip == n
            skip = n - j
        end
        if bloom_bits <= 32 # argument: is near max(p(bloom_skip)*bloom_skip)
            hash = _search_bloom_mask(tj)
            if hash&bloom_mask == 0
                if bloom_bits < 32
                    # put in bloom filter up to 32 bits
                    bloom_mask |= hash
                end
                bloom_bits += 1 # gets >32 to terminate bloom filter filling
            end
            if bloom_bits <= 32
                bloom_skip += 1
            end
        else
            if skip<n
                # we have finished adding hashes to bloom filter
                # and we have determined the skip distance if matching last byte
                # nothimg remains to be done in preprocessing - stop work.
                break
            end
        end
    end
    return t,bloom_mask,bloom_skip,bloom_bits,skip,tlast
end

function bloom_besti(s::SearchSequence, t::SearchSequence, i::Integer,sv::MaybeVector=nothing)
    bloom_besti(s,bloom_besti(t),i,sv)
end


function bloom_besti(s::SearchSequence,p::Tuple,i::Integer,sv::MaybeVector=nothing)
    @inbounds begin
    (t,bloom_mask,bloom_skip,bloom_bits,skip,tlast) = p
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

    tlast = _nthbyte(t,n)
    if DOSTATS loops = 0 end
    if DOSTATS bloomtests = 0 end
    if DOSTATS bloomskips = 0 end
    i +=n-1
    if bloom_bits <= 12 # best guess from benchmarks
        # do bloom test first in loop
        while i < m
            if DOSTATS loops += 1 end
            if DOSTATS bloomtests += 1 end
            if bloom_mask & _search_bloom_mask(_nthbyte(s,i)) == 0
                if DOSTATS bloomskips += 1 end
                i += bloom_skip
            elseif _nthbyte(s,i) == tlast
                # check candidate
                j = 1
                while j < n
                    if _nthbyte(s,i-n+j) != _nthbyte(t,j)
                        break
                    end
                    j += 1
                    # match found?
                    if j == n
                        if DOSTATS sv[Int(SFloops)] = loops; sv[Int(SFtests)] = bloomtests; sv[Int(SFskips)] = bloomskips; sv[Int(SFbits)] = bitcount(bloom_mask) end
                        return i+1
                    end
                end
                # no match: skip and test bloom
                i += skip
            else
                i +=1
            end
        end
    else
        # do byte test first
        while i < m
            if DOSTATS loops += 1 end
            if _nthbyte(s,i) == tlast
                # check candidate
                j = 1
                while j < n
                    if _nthbyte(s,i-n+j) != _nthbyte(t,j)
                        break
                    end
                    j += 1
                    # match found?
                    if j == n
                        if DOSTATS sv[Int(SFloops)] = loops; sv[Int(SFtests)] = bloomtests; sv[Int(SFskips)] = bloomskips; sv[Int(SFbits)] = bitcount(bloom_mask) end
                        return i+1
                    end
                end
                # no match: skip and test bloom
                i += skip
                if DOSTATS bloomtests += 1 end
                if i<w && bloom_mask & _search_bloom_mask(_nthbyte(s,i)) == 0
                    if DOSTATS bloomskips += 1 end
                    i += bloom_skip
                end
            else
                i += 1
                if DOSTATS bloomtests += 1 end
                if bloom_mask & _search_bloom_mask(_nthbyte(s,i)) == 0
                    if DOSTATS bloomskips += 1 end
                    i += bloom_skip
                end
            end
        end
    end
    if i==m
        # test end match
        j = 1
        while j <= n
            if _nthbyte(s,i-n+j) != _nthbyte(t,j)
                break # not found
            end
            if j == n
                if DOSTATS sv[Int(SFloops)] = loops; sv[Int(SFtests)] = bloomtests; sv[Int(SFskips)] = bloomskips; sv[Int(SFbits)] = bitcount(bloom_mask) end
                return i-n+1
            end # match at the very end
            j += 1
        end
    end
    if DOSTATS sv[Int(SFloops)] = loops end
    if DOSTATS sv[Int(SFtests)] = bloomtests end
    if DOSTATS sv[Int(SFskips)] = bloomskips end
    if DOSTATS sv[Int(SFbits)] = bitcount(bloom_mask) end
    0
end
end