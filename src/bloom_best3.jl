"""
bloom_best variant

    * use "bloom first" (bloom_v4) for small patterns

    * use struct for initialitation results directlye

benchmark result is: seems like optimizer cannot remove struct access,
so runtime is worse.
"""
function bloom_best3 end

struct SearchStructure
    t::SearchSequence
    bloom_mask::UInt64
    bloom_skip::Int
    bloom_bits::Int
    skip::Int
    tlast::UInt8
end


function bloom_best3(t::SearchSequence)
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
    skip -=1
    return SearchStructure(t,bloom_mask,bloom_skip,bloom_bits,skip,tlast)
end

function bloom_best3(s::SearchSequence, t::SearchSequence, i::Integer)
    bloom_best3(s,bloom_best3(t),i)
end


function bloom_best3(s::SearchSequence,p::SearchStructure,i::Integer,sv::MaybeVector=nothing)
    DOSTATS = !(sv isa Nothing)
    n = sizeof(p.t)
    m = sizeof(s)
    if n == 0
        return 1 <= i <= m+1 ? max(1, i) : 0
    elseif m == 0
        return 0
    elseif n == 1
        return something(findnext(isequal(_nthbyte(p.t,1)), s, i), 0)
    end
    w = m - n
    if w < 0 || i - 1 > w
        return 0
    end


    if DOSTATS loops = 0 end
    if DOSTATS bloomtests = 0 end
    if DOSTATS bloomskips = 0 end
    i +=n-1
    if p.bloom_bits <= 12 # best guess from benchmarks
        # do bloom test first in loop
        while i < m
            if DOSTATS loops += 1 end
            if DOSTATS bloomtests += 1 end
            if p.bloom_mask & _search_bloom_mask(_nthbyte(s,i)) == 0
                if DOSTATS bloomskips += 1 end
                i += p.bloom_skip
                continue
            end
            if _nthbyte(s,i) == p.tlast
                # check candidate
                j = 1
                while j < n
                    if _nthbyte(s,i-n+j) != _nthbyte(p.t,j)
                        break
                    end
                    j += 1
                    # match found?
                    if j == n
                        if DOSTATS recordcase(sv, loops, bloomtests, bloomskips, bitcount(bloom_mask), skip) end
                        return i-n+1
                    end
                end
                # no match: skip and test bloom
                i += p.skip
            end
            i +=1
        end
    else
        # do byte test first
        while i < m
            if DOSTATS loops += 1 end
            if _nthbyte(s,i) == p.tlast
                # check candidate
                j = 1
                while j < n
                    if _nthbyte(s,i-n+j) != _nthbyte(p.t,j)
                        break
                    end
                    j += 1
                    # match found?
                    if j == n
                        if DOSTATS recordcase(sv, loops, bloomtests, bloomskips, bitcount(bloom_mask), skip) end
                        return i-n+1
                    end
                end
                # no match: skip and test bloom
                i += p.skip
                if i>=w
                    break
                end
            end
            i += 1
            if DOSTATS bloomtests += 1 end
            if p.bloom_mask & _search_bloom_mask(_nthbyte(s,i)) == 0
                if DOSTATS bloomskips += 1 end
                i += p.bloom_skip
            end
        end
    end
    if i==m
        # test end match
        j = 1
        while j <= n
            if _nthbyte(s,i-n+j) != _nthbyte(p.t,j)
                break # not found
            end
            if j == n
                if DOSTATS recordcase(sv, loops, bloomtests, bloomskips, bitcount(bloom_mask), skip) end
                return i-n+1
            end # match at the very end
            j += 1
        end
    end
    if DOSTATS recordcase(sv, loops, bloomtests, bloomskips, bitcount(bloom_mask), skip) end
    0
end
