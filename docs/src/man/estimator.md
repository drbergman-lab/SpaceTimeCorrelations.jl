# The estimator

`C(Δx, Δy, Δt)` is the mean product of z-scored values over all agent pairs whose displacement
falls in a bin, pooled over every pair of frames `Δt` apart.

!!! tierbrief
    Values are z-scored within each frame, every agent in frame `t` is paired with every agent in
    frame `t + Δt` near it, and the products are averaged by displacement bin.

!!! tierfull
    For each lag `Δt` (in frames) and each frame `t` with `t + Δt` in range:

    1. Within each frame, values are z-scored: `zᵢ = (vᵢ - mean) / sd`.
    2. Each agent `i` in frame `t` gets a base point `bᵢ` (its own position at `t` in the lab
       frame; see [Lab frame and tissue frame](frames.md) for the alternative).
    3. Each agent `j` in frame `t + Δt` within reach of `bᵢ` contributes `zᵢ(t) · zⱼ(t + Δt)` to
       the bin containing `d = xⱼ(t + Δt) - bᵢ`.

    `C` in a bin is the summed product divided by the pair count. At lag 0 every agent pairs with
    itself in the zero bin, so for an uncorrelated value `C` is 1 there (diluted by any other
    agents in the same bin) and near 0 everywhere else.

```julia
r = spaceTimeCorrelation(positions, values; lags = 0:10, dr = 10.0, rmax = 100.0)
r.C[ix, iy, k]      # correlation at (r.centers[ix], r.centers[iy]), lag r.lags[k]
r.counts[ix, iy, k] # pairs behind it
```

!!! tierjournal "2026-10-08 — Per-frame z-scoring"
    Decided: z-score each frame separately, with the population standard deviation.
    Why: a frame-wide shift in expression (the whole sheet brightening over a run) is a trend, not
    a spatial correlation, and pooled normalisation would read it as long-range positive
    correlation at every lag. The population SD makes the self-pair term average exactly 1, so the
    zero bin at lag 0 reads 1 for uncorrelated values.
    Open: whether pooled normalisation should be an option, for cases where the frame-to-frame
    change in level is itself the signal.

## Direction and sign

The displacement points from the earlier agent to the later one, so a pattern moving in `+x`
puts the peak at positive `Δx`.

!!! tierfull
    The map is not symmetrised. At lag 0 it is symmetric anyway, since the pair `(i, j)` and the
    pair `(j, i)` land at `d` and `-d` with the same product. At later lags the asymmetry is the
    signal: the peak's offset from zero, divided by the lag, is the pattern's velocity.

!!! tierjournal "2026-10-08 — Directed displacement, not symmetrised"
    Decided: keep `d = x_j(t+Δt) - b_i` with its sign, and return the full plane.
    Rejected: folding `d` and `-d` together, which halves the noise but erases the direction of
    travel, the main thing the lagged maps are for.

## The grid

Bins are squares of width `dr`, one of them centred exactly on zero displacement, covering every
displacement with each coordinate within `rmax`.

!!! tierbrief
    `rmax` is rounded out to a whole number of bins: `centers` runs from `-m·dr` to `m·dr` with
    `m = ceil(rmax / dr)`.

!!! tierfull
    The grid is the full square, corners included. The neighbour search reaches the corners of
    the outermost bins, and pairs falling outside the square are dropped, so every bin is sampled
    the same way.

    Bins with fewer than `mincount` pairs (default 1, i.e. only empty bins) are `NaN`. Raise
    `mincount`, or mask with `r.counts`, to hide bins whose value rests on a handful of pairs.

```@example grid
using SpaceTimeCorrelations, Random
data = noiseFrames(; n = 400, nframes = 3, rng = Xoshiro(1))
r = spaceTimeCorrelation(data.positions, data.values; lags = 0:1, dr = 1, rmax = 2.5)
r.centers
```

```@example grid
r.counts[:, :, 2]   # pairs per bin at lag 1: corners are populated
```

!!! tierdev
    The search radius is `(m + ½)·dr·√dims`, the circumradius of the grid. A pair is kept when
    every coordinate of `d` is strictly inside `(m + ½)·dr`, then binned with `round(d / dr)`.
    Both tests are symmetric under `d → -d`, which is what makes the lag-0 map exactly symmetric
    in its counts; the `symmetry at lag 0` test relies on it.

!!! tierjournal "2026-10-08 — Square grid, searched to its corners"
    Decided: search the circumscribing radius and keep the square grid whole.
    Rejected: searching a disk of radius `rmax`, which leaves the corner bins permanently empty
    and makes "empty" mean two different things; and defining the output as the disk, which
    needs a mask to explain every plot. The extra search area costs a factor of about 4/π in
    pairs considered.

## What the numbers support

`C` is a descriptive statistic. Use `counts` to judge which bins to trust; there are no standard
errors.

!!! tierfull
    **Edge effects.** Near the domain boundary, large displacements have fewer partners, so the
    outer bins rest on fewer pairs. No edge correction is applied; `counts` shows where this
    bites.

    **No naive standard errors.** Pairs are not independent: each agent appears in many pairs, and
    neighbouring agents carry similar values. A standard error from the pair count would be far
    too small. Resample whole simulations (replicates) instead.

    **Uniform sampling.** Lags are in frames, so frames must be equally spaced in time.

    **Degenerate frames.** A frame with fewer than two agents, or with every value equal, has no
    z-score. It is left out of every pair it would join, and listed in `r.skipped`.

## Why not the pair correlation function or Ripley's K

!!! tierfull
    The pair correlation function `g(r)` and Ripley's K are the standard tools for spatial
    clustering, and they answer a different question: how the *points themselves* are arranged,
    or how points of one category sit relative to another. They take positions, optionally with
    labels, not a graded value.

    To use them on expression you would have to threshold it into "high" and "low" agents.
    That discards the graded signal: an agent at the 51st percentile and one at the 99th become
    the same, the answer depends on an arbitrary cutoff, and in a sheet of roughly uniformly
    packed agents the point pattern itself carries little information. `C` uses the value
    directly, so it measures how expression covaries with displacement, not how labelled points
    cluster.

!!! tierjournal "2026-10-08 — No thresholding of values"
    Decided: correlate the continuous value; never threshold it into categories.
    Rejected: thresholding expression and running a pair correlation function or Ripley's K on
    the "high" agents. The cutoff is arbitrary, the result moves with it, and the graded signal
    that distinguishes a weak patch from a strong one is thrown away.
