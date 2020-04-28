"""
brute force search without skips>1
"""
function naivesearch end

function naivesearch(t::SearchSequence)
    return (t,)
end

function naivesearch(s::SearchSequence, t::SearchSequence, i::Integer,sv::MaybeVector=nothing)
    naivesearch(s,naive(t),i,sv)
end

function naivesearch(s::SearchSequence,p::Tuple, i::Integer,sv::MaybeVector=nothing)
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
    (t,) = p

    tlast = _nthbyte(t,n)
    if DOSTATS loops = 0 end
    if DOSTATS bloomtests = 0 end
    if DOSTATS bloomskips = 0 end
    i -= 1
    while i <= w
        if DOSTATS loops += 1 end
        if _nthbyte(s,i+n) == tlast
            # check candidate
            j = 1
            while j < n
                if _nthbyte(s,i+j) != _nthbyte(t,j)
                    break
                end
                j += 1
                # match found?
                if j == n
                    if DOSTATS sv[Int(SFloops)] = loops; sv[Int(SFtests)] = bloomtests; sv[Int(SFskips)] = bloomskips; sv[Int(SFbits)] = 0 end
                    return i+1
                end
            end
        end
        i += 1
    end
    if DOSTATS sv[Int(SFloops)] = loops end
    if DOSTATS sv[Int(SFtests)] = bloomtests end
    if DOSTATS sv[Int(SFskips)] = bloomskips end
    if DOSTATS sv[Int(SFbits)] = 0 end
    0
end
