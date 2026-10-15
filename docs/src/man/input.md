# Input formats

Pass either per-frame arrays of positions and values, or a vector of per-frame tables.

## Arrays

!!! tierbrief
    `positions[t]` is a vector of points or a matrix with one column per agent; `values[t]` and
    `ids[t]` are vectors with one entry per agent.

```@setup input
using SpaceTimeCorrelations
include(joinpath(pkgdir(SpaceTimeCorrelations), "test", "synthetic.jl"))   # toy-data generators
```

```@example input
using SpaceTimeCorrelations, Random
data = bumpFrames(; n = 500, L = 20.0, nframes = 3, rng = Xoshiro(3))
typeof(data.positions[1])
```

```@example input
r = spaceTimeCorrelation(data.positions, data.values; ids = data.ids,
                         lags = 0:2, dr = 0.5, rmax = 2)
size(r.C)
```

!!! tierfull
    Points can be tuples, `SVector`s or plain vectors. Only the first `dims` coordinates are
    used (default 2), so 3-vectors from a 2D simulation that stores `z = 0` work unchanged. With
    `dims = 3` the grid is a cube and `C` gains a third spatial index; nothing else changes.

    Frames must be in time order and equally spaced, since `lags` count frames.

## Tables

!!! tierbrief
    Any Tables.jl table works, one per frame: `DataFrame`s, `NamedTuple`s of vectors, CSV files.
    Name the columns with `value`, `position` and (for the tissue frame) `id`.

```@example input
tables = frameTables(data)   # one NamedTuple of columns per frame
r = spaceTimeCorrelation(tables; value = :value, position = (:x, :y), id = :id,
                         lags = 0:2, dr = 0.5, rmax = 2)
r.frame
```

!!! tierfull
    `position` is either one column whose entries are points, such as PhysiCell's `position`
    column of `(x, y, z)` vectors, or a tuple of coordinate columns like `(:x, :y)`. With
    DataFrames loaded, a vector of `DataFrame`s goes in directly:

```julia
using DataFrames
frames = [DataFrame(snapshot) for snapshot in snapshots]   # one per saved time point
r = spaceTimeCorrelation(frames; value = :NICD, position = :position, id = :ID,
                         lags = 0:10, dr = 10.0, rmax = 100.0)
```

!!! tierjournal "2026-10-08 — Tables come one per frame"
    Decided: the table method takes a vector of per-frame tables.
    Rejected: a single long table with a time column. It would have to infer frame spacing from
    the time values, and an empty frame (no agents) cannot appear in it, silently shifting every
    later lag by one. Per-frame tables are also what a loop over saved snapshots produces.
