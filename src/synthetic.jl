"""
    noiseFrames(; n = 400, L = 20.0, nframes = 5, rng = Random.default_rng())

Generate fixed agents carrying independent standard-normal values: the null case, where `C` is
1 in the zero bin at lag 0 and near 0 everywhere else.

Agents are placed uniformly in the square `[0, L)²` and do not move. Returns a `NamedTuple`
`(positions, values, ids)` ready for [`spaceTimeCorrelation`](@ref).

# Examples
```jldoctest
julia> data = noiseFrames(; n = 50, nframes = 3);

julia> length(data.positions), length(data.values[1]), data.ids[1][1:3]
(3, 50, [1, 2, 3])
```
"""
function noiseFrames(; n::Integer = 400, L::Real = 20.0, nframes::Integer = 5,
                     rng::AbstractRNG = default_rng())
    x = [SVector(L * rand(rng), L * rand(rng)) for _ in 1:n]
    return (positions = [copy(x) for _ in 1:nframes],
            values = [randn(rng, n) for _ in 1:nframes],
            ids = [collect(1:n) for _ in 1:nframes])
end

"""
    bumpFrames(; n = 2000, L = 40.0, nframes = 6, σ = 1.5, density = 0.05,
               relative_velocity = (0.0, 0.0), flow = (0.0, 0.0), turnover = 0.0,
               noise = 0.0, rng = Random.default_rng())

Generate agents whose values are a field of Gaussian bumps of width `σ` that moves at
`relative_velocity` relative to the agents, while the agents themselves move with `flow`.

Every agent has a material coordinate ξ, its starting position, drawn uniformly in `[0, L)²`.
Its value at frame `t` (counting from 0) is

    v(t) = Σₘ exp(-|ξ - relative_velocity·t - cₘ|² / 2σ²) + noise · ε

with `round(density · L²)` bump centers `cₘ` placed uniformly and the field periodic with period
`L`, so the pattern has no edges. The value is therefore carried by the agent, plus whatever
motion the pattern has relative to the agents.

Positions start at ξ and advance by `flow` once per frame. `flow` is either a constant velocity
or a function of the current position returning one, which allows nonuniform motion such as two
opposing halves.

Each frame after the first, a fraction `turnover` of the agents is replaced: a newcomer with a
fresh ID takes the place, and the material coordinate, of the agent it replaces.

The ground truth, for lag `k`:

- at lag 0, the correlation at displacement `r` is `exp(-|r|² / 4σ²)`;
- the tissue-frame peak sits at `relative_velocity · k`;
- under a uniform flow `u`, the lab-frame peak sits at `(u + relative_velocity) · k`.

Returns a `NamedTuple` `(positions, values, ids)`.

# Examples
```jldoctest
julia> data = bumpFrames(; n = 500, nframes = 3, flow = (1.0, 0.0));

julia> Tuple(data.positions[3][1] - data.positions[1][1])   # two frames of flow
(2.0, 0.0)
```
"""
function bumpFrames(; n::Integer = 2000, L::Real = 40.0, nframes::Integer = 6, σ::Real = 1.5,
                    density::Real = 0.05, relative_velocity = (0.0, 0.0), flow = (0.0, 0.0),
                    turnover::Real = 0.0, noise::Real = 0.0, rng::AbstractRNG = default_rng())
    0 <= turnover <= 1 || throw(ArgumentError("turnover must be in [0, 1], got $turnover"))
    L = Float64(L)
    u = SVector{2,Float64}(relative_velocity...)
    velocity = flow isa Function ? x -> SVector{2,Float64}(flow(x)...) :
               (c = SVector{2,Float64}(flow...); _ -> c)
    centers = [SVector(L * rand(rng), L * rand(rng)) for _ in 1:round(Int, density * L^2)]
    ξ = [SVector(L * rand(rng), L * rand(rng)) for _ in 1:n]
    x = copy(ξ)
    id = collect(1:n)
    nextid = n + 1
    nreplace = round(Int, turnover * n)

    positions = Vector{Vector{SVector{2,Float64}}}(undef, nframes)
    values = Vector{Vector{Float64}}(undef, nframes)
    ids = Vector{Vector{Int}}(undef, nframes)
    for t in 0:nframes-1
        if t > 0
            x = [p + velocity(p) for p in x]
            for i in randperm(rng, n)[1:nreplace]
                id[i] = nextid
                nextid += 1
            end
        end
        shift = u * t
        values[t+1] = [_bumpfield(q - shift, centers, σ, L) + noise * randn(rng) for q in ξ]
        positions[t+1] = copy(x)
        ids[t+1] = copy(id)
    end
    return (positions = positions, values = values, ids = ids)
end

# The periodic bump field at a point, using the nearest image of each bump.
function _bumpfield(q, centers, σ, L)
    s = 0.0
    for c in centers
        d = q - c
        d = d - L * round.(d / L)
        s += exp(-(d[1]^2 + d[2]^2) / (2σ^2))
    end
    return s
end

"""
    frameTables(data) -> Vector{NamedTuple}

Turn generator output into one table per frame, with columns `id`, `x`, `y` and `value`, for
exploring the table method of [`spaceTimeCorrelation`](@ref).

Each table is a `NamedTuple` of vectors, which any Tables.jl consumer accepts; `DataFrame.(…)`
converts them.

# Examples
```jldoctest
julia> tables = frameTables(noiseFrames(; n = 10, nframes = 2));

julia> length(tables), keys(tables[1]), length(tables[1].x)
(2, (:id, :x, :y, :value), 10)
```
"""
function frameTables(data)
    return [(id = data.ids[t], x = first.(data.positions[t]), y = last.(data.positions[t]),
             value = data.values[t]) for t in eachindex(data.positions)]
end
