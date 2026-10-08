# AGENTS.md — SpaceTimeCorrelations.jl

One estimator: the joint space-time correlation `C(Δx, Δy, Δt)` of a continuous per-agent value
on moving, tracked agents (agent-based model output). The entry point is `spaceTimeCorrelation`,
with an array method in [`src/estimator.jl`](src/estimator.jl) and a per-frame-tables method in
[`src/tables.jl`](src/tables.jl) that calls it. It returns a `SpaceTimeCorrelationResult`
([`src/result.jl`](src/result.jl)). Toy-data generators with known answers live in
[`src/synthetic.jl`](src/synthetic.jl).

Start with [`docs/src/dev/architecture.md`](docs/src/dev/architecture.md): module map, the data
flow through one call, invariants. Design decisions are `!!! tierjournal` blocks on the manual
pages, collected into `docs/src/dev/journal.md`. The original brief is [`HANDOFF.md`](HANDOFF.md).

## Commands

Always run Julia with a project flag; never edit `Manifest.toml` or add a dependency without asking.

```bash
julia --project=. -e 'using Pkg; Pkg.test()'                  # the test suite
julia --project=. -t 4 -e 'using Pkg; Pkg.test()'             # same, exercising the threads
julia --project=docs -e 'using Pkg; Pkg.instantiate()'         # once; docs/Project.toml sources ".."
julia --project=docs docs/make.jl                             # build docs into docs/build, runs doctests
julia --project=docs -e 'include("docs/journal.jl")'          # regenerate docs/src/dev/journal.md alone
```

The docs build needs at least one commit (Documenter reads `HEAD`), not a remote.

## Invariants

- **Scope is one estimator.** No plotting, no significance tests, no PCF/Ripley's K, no
  wavelets/DMD/optical flow, no PhysiCell-specific I/O. Plotting dependencies are out.
- **Runtime dependencies stay minimal:** NearestNeighbors, Tables, Statistics, StaticArrays,
  Random. GeoStatsFunctions, GeoTables and Meshes are test-only (`[extras]`) and must stay so.
- **Directed displacement** `d = x_j(t + Δt) - b_i`; never symmetrise.
- **IDs locate the base agent only.** In the tissue frame `b_i` is agent i's own position at
  `t + Δt`; partners are found by position in frame `t + Δt`. Never replace this with a global
  shift.
- **Zero is a bin centre; the grid is square and searched to its corners.**
- **Per-frame z-score with the population SD.** Degenerate frames are skipped and reported in
  `skipped`, never silently dropped.
- **`lags`, `dr`, `rmax` are required.** Do not invent defaults for scales.
- The inner loop (`_accumulate!`) must not allocate. Per-task accumulators come from
  `_accumulatechunk`, a separate function: a closure that assigns `S`/`N` would capture the outer
  variables and silently share them between tasks.

## Conventions

- Functions `camelCase` (`spaceTimeCorrelation`), internal helpers `_camelCase` or `_lowercase`,
  types `PascalCase`, files `snake_case.jl`.
- Every exported name has a docstring with a runnable `jldoctest`. `checkdocs = :all`.
- Docs follow the tiered-docs house style: `!!! tierbrief` / `tierfull` / `tierdev` /
  `tierjournal`, shallowest first within a section, never wrapping a code block.
  `docs/src/dev/journal.md` is generated; edit the blocks on the pages.

## Renaming the package

The name is a placeholder. To rename: change `name` in `Project.toml`, rename
`src/SpaceTimeCorrelations.jl` and the `module` line, then replace the string everywhere:

```bash
grep -rl SpaceTimeCorrelations --exclude-dir=.git --exclude-dir=build . | xargs sed -i '' 's/SpaceTimeCorrelations/NewName/g'
```

Keep the UUID.

## Sharp edges

- `noiseFrames`/`bumpFrames` turnover keeps the replaced agent's position and value and changes
  only the ID. That is deliberate: it isolates the effect of losing base agents.
- `ntasks` changes the summation order, so `C` can differ in the last bits between task counts;
  `counts` never does.
- Lab and tissue frames coincide exactly when agents do not move (test 3 asserts `isequal`).
