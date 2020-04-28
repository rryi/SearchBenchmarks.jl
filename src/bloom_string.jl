"""
variant: test bloomfilter first
"""
function _searchindex_string(s::String, t::String, i::Integer,sv::MaybeVector==nothing)
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

    t_ptr = pointer(t)
    t_last = t_ptr + n - 1
    s_ptr = pointer(s)
    s_last = s_ptr + m - 1
    GC.@preserve s t begin
        skip = n
        tlast = unsafe_load(t_last)
        bloom_mask = UInt64(_search_bloom_mask(tlast))
        bloom_size = 1 # number of bytes to skip if byte not in filter
        bloom_bits = 1 # number of bits set in bloom_mask (if <= 32)
        bloom_skip = 1 # no. of bytes to skip if byte not in filter
        t_j = t_last # t_last is already processed
        while (t_j-=1)>=t_ptr
            tj = unsafe_load(t_j)
            if tj == tlast && skip == n
                skip = t_last - t_j
            end
            if bloom_bits <= 32
                # 32 bits in bloom filter is guess of max(expected skip per test)
                # simplified argumentation:
                # assume no hash collisions for pattern bytes.
                # then we have bloom_bits==bloom_skip.
                # assume uniform distribution of bytes to test. Then,
                # p(byte in filter) == (64-bloom_bits)/64
                # Expected skip is p(byte in filter)*bloom_skip
                # writing b for bloom_skip(==bloom_bits) we have
                # expectedSkip(b)= (64-b)/64*b
                # expectedSkip has its maximum at 32.
                hash = _search_bloom_mask(tj)
                if hash&bloom_mask == 0
                    if bloom_bits < 32
                        # put in bloom filter up to 32 bits
                        bloom_mask |= hash
                    end
                    bloom_bits += 1 # gets >32 to terminate bloom filter filling
                end
                if bloom_bits <= 32
                    # enlarge bloom_skip as long as all hashes are added to filter.
                    # while bloom_bits<=32, all hashes were added to filter
                    # or a hash collision happened (adding hash does not change filter)
                    bloom_skip += 1
                else
                    if skip<n
                        # we have finished adding hashes to bloom filter
                        # and we have determined the skip distance if matching last byte
                        # nothimg remains to be done in preprocessing - stop work.
                        break
                    end
                end
            end
        end

        if DOSTATS loops = 0 end
        if DOSTATS bloomtests = 0 end
        if DOSTATS bloomskips = 0 end
        si = s_ptr+n-1
        if bloom_bits <= 12 # best guess from benchmarks
            # do bloom test first in loop
            while si < s_last
                if DOSTATS loops += 1 end
                if DOSTATS bloomtests += 1 end
                silast = unsafe_load(si)
                if bloom_mask & _search_bloom_mask(silast) == 0
                    if DOSTATS bloomskips += 1 end
                    si += bloom_skip
                elseif silast == tlast
                    # check candidate
                    j = n
                    while (j-=1) > 0 # case j==0 is already tested (silast==tlast)
                        if unsafe_load(si-j) != unsafe_load(t_last-j)
                            break
                        end
                        # match found?
                        if j == 1
                            if DOSTATS sv[Int(SFloops)] = loops; sv[Int(SFtests)] = bloomtests; sv[Int(SFskips)] = bloomskips; sv[Int(SFbits)] = bitcount(bloom_mask) end
                            return si-s_ptr-n+2 # si points to last byte of pattern
                        end
                    end
                    # no match: skip and test bloom
                    si += skip
                else
                    si +=1
                end
            end
        else
            # do byte test first
            while si < s_last
                if DOSTATS loops += 1 end
                if  unsafe_load(si) == tlast
                    # check candidate
                    j = n
                    while (j-=1) > 0
                        if unsafe_load(si-j) != unsafe_load(t_last-j)
                            break
                        end
                        # match found?
                        if j == 1
                            if DOSTATS sv[Int(SFloops)] = loops; sv[Int(SFtests)] = bloomtests; sv[Int(SFskips)] = bloomskips; sv[Int(SFbits)] = bitcount(bloom_mask) end
                            return si-s_ptr-n+2
                        end
                    end
                    # no match: skip and test bloom
                    si += skip
                    if DOSTATS bloomtests += 1 end
                    if si<s_last && bloom_mask & _search_bloom_mask(unsafe_load(si)) == 0
                        if DOSTATS bloomskips += 1 end
                        si += bloom_skip
                    end
                else
                    # most frequent case
                    si += 1
                    if DOSTATS bloomtests += 1 end
                    if bloom_mask & _search_bloom_mask(unsafe_load(si)) == 0
                        if DOSTATS bloomskips += 1 end
                        si += bloom_skip
                    end
                end
            end
        end
        # main loop may end before last match candidate
        if si==s_last
            # test end match
            j = n
            while (j-=1) >= 0
                if unsafe_load(si-j) != unsafe_load(t_last-j)
                    break
                end
                # match found?
                if j == 0
                    if DOSTATS sv[Int(SFloops)] = loops; sv[Int(SFtests)] = bloomtests; sv[Int(SFskips)] = bloomskips; sv[Int(SFbits)] = bitcount(bloom_mask) end
                    return si-s_ptr-n+2
                end
            end
        end
    end # GC.preserve
    if DOSTATS sv[Int(SFloops)] = loops end
    if DOSTATS sv[Int(SFtests)] = bloomtests end
    if DOSTATS sv[Int(SFskips)] = bloomskips end
    if DOSTATS sv[Int(SFbits)] = bitcount(bloom_mask) end
    0
end
