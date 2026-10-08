# SpaceTimeCorrelations.jl

Measure how a continuous per-agent value is correlated across space and time in agent-based model
output, resolved by direction: `C(Δx, Δy, Δt)`.

!!! tierbrief
    Give it per-frame positions and values (and IDs, if agents are tracked); it returns a grid of
    correlations, one 2D map per time lag.

!!! tierfull
    The motivating case is NICD expression in a stem-cell sheet: is high expression patchy, how
    large are the patches, and do the patches move? A single map answers the first two at lag 0.
    Following the peak of the map across lags answers the third, and in which direction.

    Agents move, so "the patches move" has two meanings. In the **lab frame** a drift includes the
    agents' own motion; in the **tissue frame** (available when agents carry IDs) it is the motion
    of the pattern through the agents. See [Lab frame and tissue frame](man/frames.md).

```@example index
using SpaceTimeCorrelations, Random

# toy data: Gaussian bumps of expression drifting +1 per frame in x through fixed agents
data = bumpFrames(; n = 2000, L = 30.0, nframes = 4, relative_velocity = (1.0, 0.0),
                  rng = Xoshiro(1))

r = spaceTimeCorrelation(data.positions, data.values; ids = data.ids,
                         lags = 0:3, dr = 0.5, rmax = 4)
```

## Reading the result

The peak moves one unit per frame along `Δx`, the drift built into the toy data.

!!! tierbrief
    `r.C[ix, iy, k]` is the correlation at displacement `(r.centers[ix], r.centers[iy])` and lag
    `r.lags[k]`; `r.counts` holds the number of pairs behind each value.

```@example index
m = (length(r.centers) + 1) ÷ 2               # the bin centred on zero displacement
[r.centers[argmax(r.C[:, m, k])] for k in eachindex(r.lags)]   # peak along Δx at each lag
```

## Where to go next

- [The estimator](man/estimator.md): what is computed, the grid, and what the numbers mean.
- [Lab frame and tissue frame](man/frames.md): what IDs change.
- [Input formats](man/input.md): arrays, matrices, DataFrames and other tables.
- [Toy data](man/toydata.md): generators with known answers, for exploring the method.
- [API reference](reference.md).
