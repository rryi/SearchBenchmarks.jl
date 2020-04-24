using Documenter, SearchBenchmarks

makedocs(;
    modules=[SearchBenchmarks],
    format=Documenter.HTML(),
    pages=[
        "Home" => "index.md",
    ],
    repo="https://github.com/rryi/SearchBenchmarks.jl/blob/{commit}{path}#L{line}",
    sitename="SearchBenchmarks.jl",
    authors="Robert Rudolph",
    assets=String[],
)

deploydocs(;
    repo="github.com/rryi/SearchBenchmarks.jl",
)
