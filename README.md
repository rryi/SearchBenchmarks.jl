# SearchBenchmarks

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://rryi.github.io/SearchBenchmarks.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://rryi.github.io/SearchBenchmarks.jl/dev)
[![Build Status](https://travis-ci.com/rryi/SearchBenchmarks.jl.svg?branch=master)](https://travis-ci.com/rryi/SearchBenchmarks.jl)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/rryi/SearchBenchmarks.jl?svg=true)](https://ci.appveyor.com/project/rryi/SearchBenchmarks-jl)
[![Codecov](https://codecov.io/gh/rryi/SearchBenchmarks.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/rryi/SearchBenchmarks.jl)
[![Build Status](https://api.cirrus-ci.com/github/rryi/SearchBenchmarks.jl.svg)](https://cirrus-ci.com/github/rryi/SearchBenchmarks.jl)

## Status

Package compiles passes tests and builds at least with Windows 10 and Julia 1.4.1. I have no idea why CI tools report "unable to find commit in GitHub" or "build failed". If someone has an idea - please help. I would really like to get CI tools running as expected, but I am new to GitHub and linked CI tools  - this is my first Julia package in a "ready to use" state (I would not claim release 1.0, though).

## Purpose

SearchBenchmarks collects a couple of search algorithms and evaluates them in a simple benchmark framework. It compares the Julia base ByteVector search with several variants and other search algorithms, currently naive search and a Boyer Moore variant. The final goal is to discuss and recommend a new/modified implementation for the search in Julia base.

You are kindly invited to participate. I will add all search functions to the package which are sent to me and pass package tests.

## Installation


On the Pkg command line:

    (v1.0) pkg> add Random
    (v1.0) pkg> add Plots
    (v1.0) pkg> develop "https://github.com/rryi/SearchBenchmarks.jl.git"

## Usage

Benchmarking is done like this:

    using SearchBenchmarks
    benchmark(100, 1000000, "0x7f,1234", [naivesearch,bloom_v0], "julia search vs naive search on ASCII alphabet with random data,false)

This will generate a randon sequence of 1000000 bytes with seed 1234 and search with patterns at the end of the search sequence of length up to 100. See also benchmark doc string.

Instead of random data you can use files as search sequences, supply a file path as 3rd parameter. The fourth parameter is also a file path, but without extension. Path syntax is platform specific, e.g. under Windows "C:\\MyBenchmark\\Run1\\result".

However, customizing benchmarks by parameters is quite restricted. This is why I recommend to install the package in development state: start experimenting with your own benchmark settings by modifying the code. And experiment with modified or additional search functions.

## Implementing additional search functions

To work with current benchmark method and package tests, a search function must obey the following contract:

 * Function name begins with "bloom" or ends with "search" (needed for automatic include in pachage tests)

 * Method with signature (pattern::SearchSequence): returning a Tuple of preprocessed pattern properties. This method should
 contain all the initialization code before the main search loop.

 * Method with signature (searchable::SearchSequence,init::Tuple,i::Integer,sv::MaybeVector=nothing)

This method implemens the "search kernel". Results of initialization code are given in init parameter. i is the index to start search, and sv is nothing or a Vector{Int} to return statistics measures. Code segments collecting statistics must test sv for nothing. This will eliminate statistics code on pure (timing) benchmark runs.

## Current benchmark results

I have set up and run a test script analysis_?.jl in three variants, ? is TXT, BIN or RND. It compares the original julia search method **\_searchindex(s::Union{AbstractString,ByteArray}, t::Union{AbstractString,AbstractChar,Int8,UInt8}, i::Integer)** in file **base/string/search.jl** with naive search, boyer moore search and several variants which are attempts to improve runtime behaviour.

TXT is a scenario searching in redundant Utf8 text, here some julia source code. BIN searches in a binary file, here a 7zip DLL. RND uses random bytes, uniformly distributed

### Julia Base implementation - better than expected

When I was looking for an efficient search algorithm in julia, I started with a look into the Base package, expecting a lowlevel C-code implementation of a naive search (naive means here a search with comparing byte by byte, and advancing the search position by 1 on a mismatch). I was glad to see that julia took a higher order approach, and it pays off:

![graph](https://github.com/rryi/SearchBenchmarks.jl/blob/master/test/RND/(1)%20RND%20julia%20impl%20much%20better%20than%20naive%20search.png)

### small optimizations on julia code

I took a closer look on the julia implementation. I knew about bloom filters, but from the scenario of searching for the first match from a list of search phrases. Usage for simple search was new to me, and a quick web search did not reveal other cases using bloom filters for simple search. Seems to be not that common.

The first peculiarity I stumbled upon was the treatment of the "good prefix character": initialization code determines a safe skip distance in the case the last byte matches but some other byte does not, in variable _skip_. However, the implementation mostly ignores it and performs a bloom test at position i+1 instead of testing a i+skip+1. Only if bloom test fails (i.e. cannot exclude a match), skip distance is used.

My suggestion is to use the safe skip distance before applying bloom filter: I see no obvious reason why a bloom test at i+1 has better chances than a bloom test at i+skip+1. if we assume the same probability _p_ for a successful bloom test, we improve the expected skip distance by the value _(1-p)_*_skip_. However, the measurable effect is in most cases quite small, due to the low frequency of byte matches: below 0.5% on uniform distributed bytes.

The change eliminates one comparison and conditional jump, and allows for one more optimization: the two bloom tests can be combined, eleminating another comparison. Both is implemented in bloom_v1.

Next optimization addresses the _i<=w_ comparison before the bloom test: if we terminate the main loop on _i<w_ instead of _i<=w_, we can eliminate one comparison (and conditional branch) per main loop cycle. The price is an extra treatment of the case _i__w_. However in this final test we need only to compare on a complete match, no bloom test is needed. Code gets a little bit larger, but performance is improved for some percent. Implemented in bloom_v2.

Some technical improvements in the index variables are implemented in bloom_v3, removing several index operations.

The performance results of all three changes:

![graph]https://github.com/rryi/SearchBenchmarks.jl/blob/master/test/RND/(2)%20RND%20julia%20bloom%20impl%20optimizations.png)

### try bloom test first

Let-s assume an uniform byte distribution for a mowent. Then, the probability of a match of the last byte is about 0.4% (_1/256_). The probability of a successful bloom test is at least _1-n/256_ for a pattern of length n. The expected skip distance of a bloom test is then at least _n*(1-n/256)_. For n=10, we get about 9.6, which is much more than the expected skip distance of the byte test (about 1).

This considerations lead me to the idea to do the bloom test first, and the last byte test only if the bloom test does not suceed. My expectation was a clear gain, because a bloom test is quite cheap (one shift, one AND operation more than the last byte test). The idea is coded in bloom_v4.

The results were really surprising to me: outcome is highly sensitive to the pattern length and byte distribution. And the RND plot shows an unexpected local maxium for "bloom test first" at a pattern length of 32. I found no explanation for that, so far.

![graph](https://github.com/rryi/SearchBenchmarks.jl/blob/master/test/RND/(3)%20RND%20do%20bloom%20test%20first.png)

![graph](https://github.com/rryi/SearchBenchmarks.jl/blob/master/test/BIN/(3)%20BIN%20do%20bloom%20test%20first.png)

![graph](https://github.com/rryi/SearchBenchmarks.jl/blob/master/test/TXT/(3)%20TXT%20do%20bloom%20test%20first.png)

### bloom degeneration on large search patterns

With increasing pattern size, more and more bits are set in the bloom filter, reducing the probability of a successful bloom test. 
Assume again uniform byte distribution.



The question is:  Even if we take into account   Doing the last byte test,

For a pattern of length up to 10 we get about 96%. The expected skip result of the last byte test is very near to 1, the expected skip result of a bloom test is _n*(1-n/256)_  (is  in the search
On small patterns, the probability of a successful bloom test is quite high
