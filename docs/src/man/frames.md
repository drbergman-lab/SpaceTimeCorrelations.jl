# Lab frame and tissue frame

Passing agent IDs switches displacements from the fixed coordinates of the domain to each agent's
own later position, so a drift in the map is the pattern's motion through the agents rather than
including the agents' motion.

!!! tierbrief
    No `ids`: lab frame. With `ids`: tissue frame. `r.frame` records which was used.

!!! tierfull
    **Lab frame.** The base point for agent `i` is its position at `t`. A pattern that is simply
    carried along by flowing agents shows up as a drifting peak, at the agents' displacement.

    **Tissue frame.** The base point is agent `i`'s own position at `t + Δt`, found by its ID.
    The agents' motion is subtracted agent by agent, so a pattern carried by the agents keeps its
    peak at zero, and any drift that remains is motion of the pattern relative to the agents.

    Each agent carries its own displacement. Nothing assumes the flow is uniform: if two halves of
    a sheet slide in opposite directions, each base agent's own displacement is removed, and the
    tissue-frame map stays centred while the lab-frame map splits into two offset peaks.

Two halves of a sheet sliding in opposite directions, carrying their values: along `Δx` at lag
3, the lab-frame map has two peaks at `±3`, the tissue-frame map one at 0.

```@setup frames
using SpaceTimeCorrelations
include(joinpath(pkgdir(SpaceTimeCorrelations), "test", "synthetic.jl"))   # toy-data generators
```

```@example frames
using SpaceTimeCorrelations, Random
L = 30.0
flow = x -> x[2] < L / 2 ? (1.0, 0.0) : (-1.0, 0.0)   # opposing halves; values ride along
data = bumpFrames(; n = 2000, L, nframes = 4, σ = 1.0, density = 0.1, flow, rng = Xoshiro(2))

lab    = spaceTimeCorrelation(data.positions, data.values; lags = 0:3, dr = 0.5, rmax = 4)
tissue = spaceTimeCorrelation(data.positions, data.values; ids = data.ids,
                              lags = 0:3, dr = 0.5, rmax = 4)

m = (length(lab.centers) + 1) ÷ 2
k = 4   # lag 3
(lab_axis = round.(lab.C[:, m, k]; digits = 2), tissue_axis = round.(tissue.C[:, m, k]; digits = 2))
```

## IDs locate the base agent only

!!! tierfull
    The ID is used for one thing: finding where agent `i` is at `t + Δt`. Its partners `j` are
    found by position in frame `t + Δt` alone, and their earlier positions never enter. So the
    tissue frame needs only that each base agent can be found later; it does not need every
    agent tracked, and new agents can be partners as soon as they appear.

    **Turnover.** A base agent absent at `t + Δt` (it died, left, or was relabelled) contributes
    nothing at that lag. Long lags are therefore weighted toward agents that survive that long,
    and their pair counts fall. If survival correlates with the value (high-NICD cells dividing
    or dying faster, say), long-lag tissue-frame maps describe the survivors, not the sheet.

!!! tierdev
    The lookup is a `Dict` from ID to row, one per frame, built once. IDs must be unique within
    a frame (duplicates throw), and may be any hashable type. The lab and tissue frames share
    one kernel, `_accumulate!` in `src/estimator.jl`; the only branch is how `base` is chosen.

!!! tierjournal "2026-10-08 — IDs used only for the base agent"
    Decided: the tissue frame uses agent `i`'s own later position as the base point, and matches
    partners by position at `t + Δt` only.
    Rejected: estimating one global shift per frame pair (mean agent displacement) and
    subtracting it. That is exact only for uniform flow; under shear or opposing flows it leaves
    a spurious drift that looks like pattern motion. Also rejected: requiring partners to be
    tracked too, which would drop every newborn agent and adds nothing, since the partner's
    earlier position plays no role in `d`.
