# Architecture

## Module map

| File | Holds |
|:--|:--|
| `src/SpaceTimeCorrelations.jl` | Module, imports, exports, includes. |
| `src/result.jl` | `SpaceTimeCorrelationResult` and its `show`. |
| `src/estimator.jl` | The array method of `spaceTimeCorrelation` and everything it calls: input checks, z-scoring, ID maps, the threaded accumulation. |
| `src/tables.jl` | The table method: pulls columns out of per-frame tables and calls the array method. |
| `src/synthetic.jl` | `noiseFrames`, `bumpFrames`, `frameTables`. |
| `ext/SpaceTimeCorrelationsPhysiCellOutputExt.jl` | The `PhysiCellSequence` method: cell-type and dead-cell filtering, then the array method. |
| `ext/SpaceTimeCorrelationsPhysiCellModelManagerExt.jl` | The `Simulation` method: builds the sequence and forwards. |
| `test/runtests.jl` | Validation scenarios 1–8, unit checks, input forms, errors. |
| `test/crosscheck_geostats.jl` | Lag-0 comparison against GeoStatsFunctions' variogram. |
| `test/physicell.jl`, `test/fixtures/physicell_output/` | The PhysiCellOutput extension, on four real snapshots. |

## One call, start to finish

`spaceTimeCorrelation(frames; value, position, id, lags, dr, rmax)`:

1. `src/tables.jl` turns each table into a position vector, a value column and an ID column,
   and calls the array method.
2. The array method validates lengths and keywords, converts each frame's positions to
   `Vector{SVector{dims,Float64}}` (`_points`), and hands off to `_correlate`, a function barrier
   so the rest is compiled for a concrete dimension.
3. `_correlate` z-scores each frame (`_zscore`, `nothing` for a degenerate frame), builds one
   `KDTree` per usable frame and, in the tissue frame, one `Dict` from ID to row per frame.
4. It lists the work, every `(lag slot, t)` with both frames usable, and deals it round-robin to
   `ntasks` tasks. Each task (`_accumulatechunk`) owns a sum array and a count array and calls
   `_accumulate!` per frame pair.
5. `_accumulate!` picks each base point, queries the tree with a reused index buffer
   (`inrange!`), and adds products into bins. Nothing in it allocates.
6. The partial arrays are summed in task order, divided, and wrapped in the result.

## Invariants

- **Zero is a bin centre.** `centers[m + 1] == 0` exactly; binning is `round(d / dr)`.
- **The grid is square and fully searched.** The search radius circumscribes the grid; pairs
  outside the square are dropped. Shrinking the radius to `rmax` brings back empty corners.
- **Directed, unsymmetrised displacement.** `d = x_j(t + Δt) - b_i`. Do not fold `d` and `-d`.
- **IDs locate the base agent only.** Partners come from the tree of frame `t + Δt`. Do not
  replace the per-agent lookup with a global shift.
- **Population SD in the z-score**, so self-pairs average exactly 1.
- **Reproducible for fixed `ntasks`.** Counts are exact for any `ntasks`; sums can differ in the
  last bits between task counts, because floating-point addition is not associative.

## Tests

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
julia --project=. -t 4 test/runtests.jl     # directly, with threads (needs the test deps in your env)
```

Each validation scenario is its own `@testset` in `test/runtests.jl`; comment out the others, or
copy one into a REPL, to run it alone. Peaks are located with a thresholded centroid (`peak` at
the top of the file) rather than `argmax`, which wobbles by a bin on broad peaks.

## Adding to it

New toy scenarios belong in `bumpFrames` as keywords, with a row in the table on
[Toy data](../man/toydata.md) and a `@testset` that checks the stated answer. A new input form
belongs in `src/tables.jl` (or as a `_points` method) and should end in a call to the array
method, so the estimator has exactly one entry point.
