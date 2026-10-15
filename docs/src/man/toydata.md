# Toy data

Two generators produce frames with a known answer, for exploring what the maps look like before
trusting them on real output. They live with the tests, in `test/synthetic.jl`, and load with
`include(joinpath(pkgdir(SpaceTimeCorrelations), "test", "synthetic.jl"))`.

!!! tierbrief
    `noiseFrames` gives the null case. `bumpFrames` gives Gaussian bumps of expression that can
    drift relative to the agents while the agents themselves flow.

!!! tierfull
    In `bumpFrames` each agent has a material coordinate ξ (where it started), and its value is
    a periodic field of Gaussian bumps of width `σ`, evaluated at `ξ - relative_velocity·t`. So
    the value is carried by the agent, and `relative_velocity` is how fast the pattern moves
    through the agents. The agents' positions advance by `flow` each frame; `flow` can be a
    function of position for nonuniform motion.

    The known answers:

    | Scenario | Settings | Expect |
    |:--|:--|:--|
    | Noise baseline | `noiseFrames()` | 1 in the zero bin at lag 0, ≈ 0 elsewhere |
    | Cluster size | `bumpFrames(; σ)` | lag 0: `C(r) ≈ exp(-r² / 4σ²)` |
    | Translating pattern | `relative_velocity = v` | peak at `v·k`, both frames |
    | Pure advection | `flow = u` | lab peak at `u·k`; tissue peak at 0 |
    | Relative motion | `flow = u, relative_velocity = v` | lab peak at `(u + v)·k`; tissue at `v·k` |
    | Diagonal motion | `relative_velocity = (a, b)` | peak at `(a, b)·k`, off both axes |
    | Opposing flows | `flow = x -> x[2] < L/2 ? (u, 0) : (-u, 0)` | lab: peaks at `±u·k`; tissue: 0 |
    | Turnover | `turnover = f` | tissue frame runs; counts fall with lag |

    The test suite checks each row.

```@setup toy
using SpaceTimeCorrelations
include(joinpath(pkgdir(SpaceTimeCorrelations), "test", "synthetic.jl"))   # toy-data generators
```

```@example toy
using SpaceTimeCorrelations, Random
σ = 1.5
data = bumpFrames(; n = 3000, L = 40.0, nframes = 1, σ, rng = Xoshiro(4))
r = spaceTimeCorrelation(data.positions, data.values; lags = [0], dr = 0.5, rmax = 4)
m = (length(r.centers) + 1) ÷ 2
[(Δx = x, C = round(r.C[i, m, 1]; digits = 2), theory = round(exp(-x^2 / 4σ^2); digits = 2))
 for (i, x) in enumerate(r.centers) if x >= 0]
```

!!! tierfull
    **Why `exp(-r² / 4σ²)`?** Bumps placed at random are a shot-noise field; its covariance is the
    overlap of two Gaussians of width `σ` a distance `r` apart, which is a Gaussian of width
    `√2·σ`. The decay is wider than a single bump, which is worth remembering when reading a
    cluster size off a real map.

    **Turnover**, as generated here, gives a newcomer the place of the agent it replaces. The
    values do not change; only the IDs do, which isolates the effect of losing base agents.
