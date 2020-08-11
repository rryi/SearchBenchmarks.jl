# SearchBenchmarks

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://rryi.github.io/SearchBenchmarks.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://rryi.github.io/SearchBenchmarks.jl/dev)
[![Build Status](https://travis-ci.com/rryi/SearchBenchmarks.jl.svg?branch=master)](https://travis-ci.com/rryi/SearchBenchmarks.jl)
[![Codecov](https://codecov.io/gh/rryi/SearchBenchmarks.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/rryi/SearchBenchmarks.jl)


## Purpose

SearchBenchmarks chooses some search algorithms and evaluates them in a simple benchmark framework. Technically, we talk about search for a byte sequence, but in julia, it also covers searching for substrings in strings. String search is widely used, its performance can have substantial impact on application performance.

SearchBenchmarks compares the runtime performance of Julia base ByteVector search with several variants and other search algorithms, currently naive search and a Boyer Moore variant. Its objective is to discuss and finally recommend a new/modified implementation for byte sequence (including string) search in Julia base.

You are kindly invited to participate. I will add all contributed search functions to SearchBenchmarks which pass my package tests.

## Current State

### technical

SearchBenchmarks compiles, passes tests and builds - at least with Windows 10 and Julia 1.4.1, 1.4.2 and 1.5.0. I have no idea why some CI tools report "unable to find commit in GitHub" or "build failed". If someone has an idea - please help. I would really like to get CI tools running as expected, but I am new to GitHub and linked CI tools  - this is my first Julia package in a "ready to use" state (I would not claim release 1.0, though).

### findings and recommendation

I suggest a couple of code changes to the current julia base implementation, improving its runtime performance in almost all tested scenarios by 20% and more in the implementation bloom_best. For heavy searching, a search function based on some Boyer-Moore algorithm seems superior, with gains of up to 50% in the search kernel, but at higher initialization costs. My current idea is to add it as an additional search function, split into search initialization and search kernel.

For details, see benchmark results.

## Installation

On the Pkg command line:

    (v1.0) pkg> add Random
    (v1.0) pkg> add Plots
    (v1.0) pkg> develop "https://github.com/rryi/SearchBenchmarks.jl.git"

## Usage

Benchmarking is done like this:

    using SearchBenchmarks
    benchmark(100, 1000000, "0x7f,1234", [naivesearch,bloom_v0], "julia search vs naive search on ASCII alphabet with random data",false)

This will generate a random sequence of 1000000 bytes in the range 0:0x7f with seed 1234 and search with patterns at the end of the search sequence of length up to 100. See also benchmark doc string.

Instead of random data, you can use files as sequences to search in - supply a String with a file path as 3rd parameter. The fourth parameter is also a file path, but without extension. Path syntax is platform specific, e.g. in Windows "C:\\MyBenchmark\\Run1\\result". The last parameter causes collection of statistics, if set to **true**. Collecting statistics distorts runtime measurement, it is intended to be used in separate runs for debugging and explaining runtime behavior.

However, customizing benchmarks by parameters is quite restricted. This is why I recommend to install the package in development state: start experimenting with your own benchmark settings by modifying the code. And experiment with modified or additional search functions.

## Implementing additional search functions

To work with current benchmark method and package tests, a search function must obey the following contract:

 * Function name begins with "bloom" or ends with "search" (needed for automatic inclusion in package tests)

 * Method with signature (pattern::SearchSequence): returning a Tuple of preprocessed pattern properties. This method should
 contain all the initialization code before the main search loop.

 * Method with signature (searchable::SearchSequence,init::Tuple,i::Integer,sv::MaybeVector=nothing)

This method implements the "search kernel". Results of initialization code are supplied with init parameter. i is the index to start search, and sv is nothing or a Vector{Int} to return statistics measures. Code segments collecting statistics must test sv for nothing. This will eliminate statistics code on pure (timing) benchmark runs.

## Current benchmark results

I have set up and run a test script analysis_?.jl in three variants, ? is TXT, BIN or RND. It compares the original julia search method **\_searchindex(s::Union{AbstractString,ByteArray}, t::Union{AbstractString,AbstractChar,Int8,UInt8}, i::Integer)** in file **base/string/search.jl** with naive search, boyer moore search and several variants which are attempts to improve runtime behaviour.

TXT is a scenario searching in redundant Utf8 text, my experiments use some julia source code. BIN searches in a binary file, I used a 7zip DLL. RND uses uniformly distributed random byte sequences.

All benchmark runs generate a CSV table and a graphic showing the runtime of the included functions in mSec in relation to the search pattern size in bytes.

### Julia base implementation - better than expected

When I was looking for an efficient search algorithm in julia, I started with a look into the Base package, expecting a lowlevel C-code implementation of a naive search (naive here means a search comparing byte by byte, and advancing the search position by 1 on a mismatch). I was glad to see that julia took a higher order approach, and it pays off (bloom_v0 is a copy of julia base \_searchindex with technical modifications for benchmarking):

![graph](https://github.com/rryi/SearchBenchmarks.jl/blob/master/test/RND/(1)%20RND%20julia%20impl%20much%20better%20than%20naive%20search.png)

![graph](https://github.com/rryi/SearchBenchmarks.jl/blob/master/test/BIN/(1)%20BIN%20julia%20impl%20much%20better%20than%20naive%20search.png)

![graph](https://github.com/rryi/SearchBenchmarks.jl/blob/master/test/TXT/(1)%20TXT%20julia%20impl%20much%20better%20than%20naive%20search.png)

### Small optimizations on julia search code

I took a closer look at the julia implementation, which is based on bloom filters. I knew about them, but in a scenario of searching for the first match from a list of search phrases. Usage for simple search was new to me, and a quick web search did not reveal other cases using bloom filters for simple search.

The first peculiarity I stumbled upon was the treatment of the "good prefix character": initialization code determines a safe skip distance in the case the last byte matches but some other byte does not, in variable _skip_. However, the implementation mostly ignores it and performs a bloom test at position i+1 instead of testing at i+skip+1. Only if bloom test fails (i.e. cannot exclude a match), skip distance is used.

My suggestion is to use the safe skip distance before applying bloom filter: I see no obvious reason why a bloom test at i+1 has better chances than a bloom test at i+skip+1. if we assume the same probability _p_ for a successful bloom test, we improve the expected skip distance by the value _(1-p)_*_skip_. The measurable effect will be quite small in most cases, due to the low frequency of byte matches: below 0.5% on uniform distributed bytes.

The change eliminates one comparison and conditional jump, and allows for some code optimization: the two bloom tests can be combined, eliminating one more comparison. Both is implemented in bloom_v1. Surprisingly, the change does not pay off in terms of runtime: in all tests, bloom_v1 is a little bit slower than bloom_v0. I have no idea what the reason for this is: with statistics on, one can see that bloom_v1 constantly performs less main loop cycles. Comments are welcome.

Next optimization addresses the _i<=w_ comparison before the bloom test: if we terminate the main loop on _i<w_ instead of _i<=w_, we can eliminate one comparison (and conditional branch) per main loop cycle. The price is some extra code for the case _i==w_. However in this final test we need only compare a complete match, no bloom test is needed. Code gets a little bit larger, but performance is improved for some percent. Implemented in bloom_v2.

Some technical improvements in the index variables are implemented in bloom_v3, removing several index operations.

The performance results of all three changes:

![graph](https://github.com/rryi/SearchBenchmarks.jl/blob/master/test/RND/(2)%20RND%20julia%20bloom%20impl%20optimizations.png)

![graph](https://github.com/rryi/SearchBenchmarks.jl/blob/master/test/BIN/(2)%20BIN%20julia%20bloom%20impl%20optimizations.png)

![graph](https://github.com/rryi/SearchBenchmarks.jl/blob/master/test/TXT/(2)%20TXT%20julia%20bloom%20impl%20optimizations.png)

### Try bloom test first

In the "normal" case of the kernel main loop, the last byte test fails and the following bloom test succeeds. Bloom test was introduced because it is much more efficient than naive search: for a small computational overhead (one shift, one AND operation), we get a much larger expected skip distance. Can we improve efficiency, if we start directly with the bloom test in the main loop, and do the byte compares only if bloom test fails? We will arrive at some more main loop cycles, but eliminate the last byte test in the "normal" case.

The idea is implemented in bloom_v4. Code is a little bit smaller and simpler, compared to all the variants above. The benchmark results are ... interesting. There is no clear winner. It depends on pattern size and redundancy structure: searches in text, binary code and random data perform quite differently. As a tendency, one can conclude that bloom_v4 performs worse on medium sized patterns. This is contradictory to my expectation at least for random data:

Let us assume bytes in pattern and sequence to search are stochastically independent and uniformly distributed. Then, the probability of a match of the last byte is about 0.4% (_1/256_). The probability of a successful bloom test is _1-k/64_ if k bits are set in the bloom filter. A guess of the expected skip distance of a bloom test is _n*(1-k/64)_, with n denoting the pattern length and k the number of bits set in the bloom filter. In the case of independent and uniformly distributed bytes, one can deduct the distribution of k for a given n and compute the exact value for the expected skip distance, but this is beyond the scope here. However to give an estimation, I consider the worst case for the skip distance estimate: k==n. The formula _n*(1-n/64)_ has a maximum at 32, in the range 13:51 all values are above 10, outside this range below 10. According to this, "bloom test first" should be best for medium sized search patterns, compared to "byte test first".

Benchmark results for random data were really surprising to me: the RND plot shows an unexpected local maximum for "bloom test first" at a pattern length of 32. I have found no explanation for that, so far. Comments are welcome.

![graph](https://github.com/rryi/SearchBenchmarks.jl/blob/master/test/RND/(3)%20RND%20do%20bloom%20test%20first.png)

![graph](https://github.com/rryi/SearchBenchmarks.jl/blob/master/test/BIN/(3)%20BIN%20do%20bloom%20test%20first.png)

![graph](https://github.com/rryi/SearchBenchmarks.jl/blob/master/test/TXT/(3)%20TXT%20do%20bloom%20test%20first.png)

### Bloom degeneration on large search patterns

With increasing pattern size, more and more bits are set in the bloom filter, reducing the probability of a successful bloom test. Finally, if all bits are set in the bloom mask, bloom tests get useless. Search performance degenerates and becomes worse than native search. The pattern size where degeneration starts varies with the test scenarios, but all scenarios and all bloom-based search implementations considered so far suffer from it. In the RND benchmark, one can nicely see the effect of setting the last bits in the bloom filter.

![graph](https://github.com/rryi/SearchBenchmarks.jl/blob/master/test/RND/(4)%20RND%20degeneration%20on%20large%20patterns.png)

![graph](https://github.com/rryi/SearchBenchmarks.jl/blob/master/test/BIN/(4)%20BIN%20degeneration%20on%20large%20patterns.png)

![graph](https://github.com/rryi/SearchBenchmarks.jl/blob/master/test/TXT/(4)%20TXT%20degeneration%20on%20large%20patterns.png)

Sure, in most applications searches with a pattern length above 100 are irrelevant. But: no one expects such a degeneration for a general purpose search algorithm. And it is quite easy to avoid it: simply stop setting more bits in bloom filter once it reaches some treshold. This will lower the skip distance of a successful bloom test in favor of a higher skip probability. In function bloom_best, initializing the bloom filter stops if 32 bits are set. The choice of 32 is motivated by the simple guess of the expected bloom skip distance in the preceding chapter, it is not claimed to be an optimal setting - but benchmarks show it is often near the optimum. The kernel of bloom_best is derived from bloom_v4, with some necessary technical adjustments. As one can see, bloom_best solves the degeneration problem:

![graph](https://github.com/rryi/SearchBenchmarks.jl/blob/master/test/RND/(5)%20RND%20limited%20bits%20in%20bloom%20filter.png)

![graph](https://github.com/rryi/SearchBenchmarks.jl/blob/master/test/BIN/(5)%20BIN%20limited%20bits%20in%20bloom%20filter.png)

![graph](https://github.com/rryi/SearchBenchmarks.jl/blob/master/test/TXT/(5)%20TXT%20limited%20bits%20in%20bloom%20filter.png)

### An alternative for heavy searching: Boyer-Moore

The so-called Boyer-Moore search algorithm (and variants) is probably the most widespread method for string search besides naive search. Compared to the bloom filter based approach, it has a higher memory consumption, in the implemented variant an int array of 256 elements for byte search. Initialization effort includes an array allocation and initialization, which outweighs search effort for many "small" searches. For this reason, Boyer-Moore is not recommended as general purpose search algorithm. Kernel performance tends to be superior to bloom based search, at least for larger search pattern sizes:

![graph](https://github.com/rryi/SearchBenchmarks.jl/blob/master/test/RND/(6)%20RND%20boyer%20moore%20alternative.png)
![graph](https://github.com/rryi/SearchBenchmarks.jl/blob/master/test/BIN/(6)%20BIN%20boyer%20moore%20alternative.png)
![graph](https://github.com/rryi/SearchBenchmarks.jl/blob/master/test/TXT/(6)%20TXT%20boyer%20moore%20alternative.png)

Boyer-Moore is not susceptible to degeneration on huge pattern sizes:

![graph](https://github.com/rryi/SearchBenchmarks.jl/blob/master/test/RND/(7)%20RND%20boyer%20moore%20on%20large%20patterns.png)
![graph](https://github.com/rryi/SearchBenchmarks.jl/blob/master/test/BIN/(7)%20BIN%20boyer%20moore%20on%20large%20patterns.png)
![graph](https://github.com/rryi/SearchBenchmarks.jl/blob/master/test/TXT/(7)%20TXT%20boyer%20moore%20on%20large%20patterns.png)

## Discussion

Starting point is the original thread when bloom based search was introduced into julia base:

<https://github.com/JuliaLang/julia/pull/3797>

Since then, no significant change was done to julia's **\_searchindex** function. I have started a new discussion here:

 <https://discourse.julialang.org/t/improvement-of-string-search-in-julia-base/41945>
