# SearchBenchmarks

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://rryi.github.io/SearchBenchmarks.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://rryi.github.io/SearchBenchmarks.jl/dev)
[![Build Status](https://travis-ci.com/rryi/SearchBenchmarks.jl.svg?branch=master)](https://travis-ci.com/rryi/SearchBenchmarks.jl)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/rryi/SearchBenchmarks.jl?svg=true)](https://ci.appveyor.com/project/rryi/SearchBenchmarks-jl)
[![Codecov](https://codecov.io/gh/rryi/SearchBenchmarks.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/rryi/SearchBenchmarks.jl)
[![Build Status](https://api.cirrus-ci.com/github/rryi/SearchBenchmarks.jl.svg)](https://cirrus-ci.com/github/rryi/SearchBenchmarks.jl)

## Warning

Package is in development.

Current version does compile, but crashes in Atom/Juno REPL. See

https://github.com/JunoLab/Juno.jl/issues/570

## Purpose

SearchBenchmarks collects a couple of search algorithms and evaluates them in a simple benchmark framework. It compares the Julia base String and ByteVector search with several variants and other search algorithms. The final goal is to discuss and recommend a new/modified implementation for the search in Julia base.

## installation


    using Pkg
    Pkg.add("Plots")
    Pkg.develop("https://github.com/rryi/SearchBenchmarks.jl.git")

or in the Pkg command line:

    (v1.0) pkg> add Plots
    (v1.0) pkg> develop "https://github.com/rryi/SearchBenchmarks.jl.git"

## usage

Package is in development - please be patient until release 1.0 is available
