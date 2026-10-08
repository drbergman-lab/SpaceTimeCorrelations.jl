using SpaceTimeCorrelations
using Documenter

# Regenerate docs/src/dev/journal.md from the `!!! tierjournal` blocks on the pages,
# and check that every page's tier blocks run shallowest first.
include("journal.jl")

DocMeta.setdocmeta!(SpaceTimeCorrelations, :DocTestSetup, :(using SpaceTimeCorrelations);
                    recursive=true)

makedocs(;
    modules=[SpaceTimeCorrelations],
    authors="Daniel Bergman <danielrbergman@gmail.com> and contributors",
    sitename="SpaceTimeCorrelations.jl",
    # Named explicitly so the docs build before the repository has a remote.
    repo=Remotes.GitHub("drbergman-lab", "SpaceTimeCorrelations.jl"),
    format=Documenter.HTML(;
        canonical="https://drbergman-lab.github.io/SpaceTimeCorrelations.jl",
        edit_link="main",
        assets=["assets/tiers.css", "assets/tiers.js"],
    ),
    checkdocs=:all,
    pages=[
        "Home" => "index.md",
        "Manual" => [
            "man/estimator.md",
            "man/frames.md",
            "man/input.md",
            "man/toydata.md",
        ],
        "API reference" => "reference.md",
        "Developers" => [
            "Architecture" => "dev/architecture.md",
            "Journal" => "dev/journal.md",
        ],
    ],
)

deploydocs(;
    repo="github.com/drbergman-lab/SpaceTimeCorrelations.jl",
    devbranch="main",
)
