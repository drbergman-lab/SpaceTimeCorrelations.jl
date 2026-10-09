using SpaceTimeCorrelations
using Documenter

# Loaded so the PhysiCellOutput extension is present and its docstring is rendered and checked.
# The PhysiCellModelManager extension is left out of the build to keep PCMM out of the docs
# environment; its docstring lives in ext/ for the REPL.
using PhysiCellOutput
const SpaceTimeCorrelationsPhysiCellOutputExt =
    Base.get_extension(SpaceTimeCorrelations, :SpaceTimeCorrelationsPhysiCellOutputExt)

# Regenerate docs/src/dev/journal.md from the `!!! tierjournal` blocks on the pages,
# and check that every page's tier blocks run shallowest first.
include("journal.jl")

DocMeta.setdocmeta!(SpaceTimeCorrelations, :DocTestSetup, :(using SpaceTimeCorrelations);
                    recursive=true)

makedocs(;
    modules=[SpaceTimeCorrelations, SpaceTimeCorrelationsPhysiCellOutputExt],
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
            "man/physicell.md",
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
