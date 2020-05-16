"""
a simplified Boyer Moore search implementation

this is basically the boyer moore horspool algorithm, which restricts itself
to the "bad character" heuristics.

If the last byte matches, the skip is taken as maximum of the
"bad char" skip and the "good char" skip like in the bloom algorithms.
From a boyer moore viewpoint, this is a simplified special case of the
"good prefix" heuristics with a 1-byte-prefix.

There are more sophisticated implementations which care for worst case
performance and a better "good character" skip.

This implementation is not meant to earn a performance crown, it shall indicate
if (and at which pattern length) bloom filter search is comparable to
boyer moore search.
"""
function bmsearch end

function bmsearch(t::SearchSequence)
    size = sizeof(t) # no preprocessing restriction on huge patterns
    skip = size # good char skip like in bloom
    tlast = _nthbyte(t,size)
    badchar = zeros(Int,256) # for small patterns, UInt8 would be faster
    for j in 1:size-1
        badchar[_nthbyte(t,j)+1] = j
        if _nthbyte(t,j)==tlast
            skip = size-j
        end
    end
    badchar[tlast+1] = size
    return (t,badchar,skip,tlast)
end

function bmsearch(s::SearchSequence, t::SearchSequence, i::Integer)
    bmsearch(s,bmsearch(t),i)
end

function bmsearch(s::SearchSequence,p::Tuple, i::Integer,sv::MaybeVector=nothing)
    (t,badchar,skip,tlast) = p
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
    if DOSTATS bloomtests = 0 end # here: badchar skip > simple goodsuffix skip
    if DOSTATS bloomskips = 0 end # here: badchar skips >1
    i += n-1
    while i <= m
        if DOSTATS loops += 1 end
        if (silast = _nthbyte(s,i)) != tlast
            i += n - badchar[silast+1] # guaranteed >0
            if DOSTATS && badchar[silast+1]>1; bloomskips += 1 end
        else
            # check candidate
            j = n-1
            while (sij=_nthbyte(s,i-n+j) == _nthbyte(t,j))
                # match found?
                if j == 1
                    if DOSTATS recordcase(sv, loops, bloomtests, bloomskips, 0, skip) end
                    return i-n+1
                end
                j -= 1
            end
            # mismatch at j
            i += max(skip,j - badchar[sij+1])
            if DOSTATS && j - badchar[sij+1]>skip; bloomtests+=1;end
        end
    end
    if DOSTATS recordcase(sv, loops, bloomtests, bloomskips, 0, skip) end
    0
end
