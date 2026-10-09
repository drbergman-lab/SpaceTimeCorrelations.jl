# SpaceTimeCorrelations.jl

Joint space-time correlation `C(Δx, Δy, Δt)` of a continuous per-agent value on moving, tracked
agents, for agent-based model output.

```julia
using SpaceTimeCorrelations

# one table per frame (DataFrames, NamedTuples, CSV.File, …)
r = spaceTimeCorrelation(frames; value = :NICD, position = :position, id = :ID,
                         lags = 0:10, dr = 10.0, rmax = 100.0)

r.C[ix, iy, k]   # correlation at displacement (r.centers[ix], r.centers[iy]), lag r.lags[k]
r.counts         # pairs behind each bin
```

With `id`, displacements are measured from each agent's own later position (tissue frame), so a
drifting peak means the pattern moves through the agents. Without it they are in fixed lab
coordinates. `noiseFrames` and `bumpFrames` generate toy data with known answers.

With `using PhysiCellOutput` (or `using PhysiCellModelManager`) it reads PhysiCell output
directly:

```julia
r = spaceTimeCorrelation(PhysiCellSequence("output"; include_cells = true);
                         value = :notch, cell_types = "rgc", lags = 0:6:60, dr = 10.0, rmax = 100.0)
```

Build the docs locally with `julia --project=docs docs/make.jl`.

## Status

v0.1, local only. The estimator, table layer, toy-data generators and validation suite are in
place. Open decisions:

- the package name (a placeholder; see `AGENTS.md` for the rename recipe);
- whether pooled normalisation becomes an option alongside per-frame z-scoring.
