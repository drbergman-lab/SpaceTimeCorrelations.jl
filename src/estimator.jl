"""
    spaceTimeCorrelation(positions, values; lags, dr, rmax, ids = nothing, dims = 2,
                         mincount = 1, ntasks = Threads.nthreads()) -> SpaceTimeCorrelationResult

Estimate the space-time correlation `C(Δx, Δy, Δt)` of a per-agent value from plain per-frame
arrays.

`positions[t]` holds the positions of the agents in frame `t`, either as a vector of points
(tuples, `SVector`s or vectors) or as a matrix with one column per agent. Only the first `dims`
coordinates are used. `values[t]` holds one value per agent. Frames are equally spaced and in
time order; `lags` are in frames.

Without `ids`, displacements are measured in the lab frame. With `ids` (`ids[t]` holds one
persistent ID per agent), they are measured in the tissue frame: from each agent's own position
at the later frame.

# Keywords
- `lags`: nonnegative frame lags, e.g. `0:10`.
- `dr`: bin width.
- `rmax`: the grid covers every displacement with each coordinate in `[-rmax, rmax]`, rounded
  out to whole bins.
- `ids`: per-frame agent IDs; selects the tissue frame.
- `dims`: number of spatial coordinates to use.
- `mincount`: bins with fewer pairs are `NaN`.
- `ntasks`: number of tasks the frame pairs are split over. The result is reproducible for a
  fixed `ntasks`.

# Examples
```jldoctest
julia> positions = [[(0.0, 0.0), (5.0, 0.0)], [(1.0, 0.0), (6.0, 0.0)]];

julia> values = [[1.0, -1.0], [1.0, -1.0]];

julia> r = spaceTimeCorrelation(positions, values; lags = [1], dr = 1, rmax = 2);

julia> r.C[:, 3, 1]   # along Δx at Δy = 0: the pattern moved +1
5-element Vector{Float64}:
 NaN
 NaN
 NaN
   1.0
 NaN
```
"""
function spaceTimeCorrelation(positions::AbstractVector, values::AbstractVector;
                              lags, dr::Real, rmax::Real, ids = nothing, dims::Integer = 2,
                              mincount::Integer = 1, ntasks::Integer = Threads.nthreads())
    nframes = length(positions)
    length(values) == nframes ||
        throw(ArgumentError("got $nframes frames of positions but $(length(values)) of values"))
    ids === nothing || length(ids) == nframes ||
        throw(ArgumentError("got $nframes frames of positions but $(length(ids)) of ids"))
    dr > 0 || throw(ArgumentError("dr must be positive, got $dr"))
    rmax > 0 || throw(ArgumentError("rmax must be positive, got $rmax"))
    dims >= 1 || throw(ArgumentError("dims must be at least 1, got $dims"))
    ntasks >= 1 || throw(ArgumentError("ntasks must be at least 1, got $ntasks"))
    lagvec = collect(Int, lags)
    all(>=(0), lagvec) || throw(ArgumentError("lags must be nonnegative, got $lags"))
    allunique(lagvec) || throw(ArgumentError("lags must be unique, got $lags"))

    pts = [_points(positions[t], Val(Int(dims)), t) for t in 1:nframes]
    for t in 1:nframes
        length(values[t]) == length(pts[t]) || throw(ArgumentError(
            "frame $t has $(length(pts[t])) positions but $(length(values[t])) values"))
        ids === nothing || length(ids[t]) == length(pts[t]) || throw(ArgumentError(
            "frame $t has $(length(pts[t])) positions but $(length(ids[t])) ids"))
    end
    return _correlate(pts, values, ids, lagvec, Float64(dr), Float64(rmax), mincount, ntasks)
end

# Positions of one frame as SVectors of the first D coordinates.
function _points(p::AbstractMatrix, ::Val{D}, t) where {D}
    size(p, 1) >= D || throw(ArgumentError(
        "frame $t: positions have $(size(p, 1)) coordinates, fewer than dims = $D"))
    out = [SVector{D,Float64}(ntuple(k -> p[k, i], Val(D))) for i in axes(p, 2)]
    return _checkfinite(out, t)
end

function _points(p::AbstractVector, ::Val{D}, t) where {D}
    out = Vector{SVector{D,Float64}}(undef, length(p))
    for (i, q) in enumerate(p)
        length(q) >= D || throw(ArgumentError(
            "frame $t: position $i has $(length(q)) coordinates, fewer than dims = $D"))
        out[i] = SVector{D,Float64}(ntuple(k -> q[k], Val(D)))
    end
    return _checkfinite(out, t)
end

function _checkfinite(out, t)
    all(q -> all(isfinite, q), out) || throw(ArgumentError("frame $t has a non-finite position"))
    return out
end

# z-score one frame with the population standard deviation, so that the mean of z² is exactly 1.
# Returns `nothing` when the z-score is undefined.
function _zscore(v, t)
    all(isfinite, v) || throw(ArgumentError("frame $t has a non-finite value"))
    length(v) < 2 && return nothing
    μ = mean(v)
    σ = std(v; mean = μ, corrected = false)
    σ > 0 || return nothing
    return Float64[(x - μ) / σ for x in v]
end

function _idmap(ids, t)
    m = Dict{eltype(ids),Int}()
    sizehint!(m, length(ids))
    for (i, id) in enumerate(ids)
        haskey(m, id) && throw(ArgumentError("frame $t: id $(repr(id)) appears more than once"))
        m[id] = i
    end
    return m
end

function _correlate(pts::Vector{Vector{SVector{D,Float64}}}, values, ids, lags, dr, rmax,
                    mincount, ntasks) where {D}
    nframes = length(pts)
    m = ceil(Int, rmax / dr)
    W = (m + 0.5) * dr          # half-width of the grid: the outer edge of the outermost bin
    radius = W * sqrt(D)        # circumscribes the grid, so corner bins are searched too
    nb = 2m + 1
    centers = collect((-m:m) .* dr)

    zs = [_zscore(values[t], t) for t in 1:nframes]
    skipped = findall(isnothing, zs)
    trees = [z === nothing ? nothing : KDTree(pts[t]) for (t, z) in enumerate(zs)]
    tissue = ids !== nothing
    idmaps = tissue ? [_idmap(ids[t], t) for t in 1:nframes] : nothing

    work = [(k, t) for k in eachindex(lags) for t in 1:nframes-lags[k]
            if zs[t] !== nothing && zs[t+lags[k]] !== nothing]
    gridsize = (ntuple(_ -> nb, Val(D))..., length(lags))
    nchunks = max(1, min(ntasks, length(work)))

    # Each task owns its accumulators. Work is dealt round-robin so that the cheap and the
    # expensive lags are spread evenly, and the partial sums are combined in chunk order.
    tasks = map(1:nchunks) do c
        Threads.@spawn _accumulatechunk(c, nchunks, work, gridsize, lags, pts, zs, trees, ids,
                                        idmaps, dr, W, m, radius)
    end
    S, N = fetch(first(tasks))
    for task in Iterators.drop(tasks, 1)
        Sc, Nc = fetch(task)
        S .+= Sc
        N .+= Nc
    end

    C = [n >= mincount && n > 0 ? s / n : NaN for (s, n) in zip(S, N)]
    return SpaceTimeCorrelationResult(C, N, centers, lags, dr, tissue ? :tissue : :lab, skipped)
end

# One task's share of the work: every nchunks-th frame pair, starting at c, into its own
# accumulators. A separate function rather than a closure, so that the accumulators are local to
# the task and the hot loop is compiled for concrete argument types.
function _accumulatechunk(c, nchunks, work, gridsize, lags, pts, zs, trees, ids, idmaps,
                          dr, W, m, radius)
    S = zeros(Float64, gridsize)
    N = zeros(Int, gridsize)
    buf = Int[]
    for w in c:nchunks:length(work)
        k, t = work[w]
        _accumulate!(S, N, buf, k, t, t + lags[k], pts, zs, trees, ids, idmaps, dr, W, m, radius)
    end
    return S, N
end

# Add every pair (i in frame t, j in frame s = t + lag) into lag slot k.
function _accumulate!(S, N, buf, k, t, s, pts, zs, trees, ids, idmaps, dr, W, m, radius)
    early, late = pts[t], pts[s]
    zt, zs_ = zs[t], zs[s]
    tree = trees[s]
    for i in eachindex(early)
        if idmaps === nothing
            base = early[i]
        else
            r = get(idmaps[s], ids[t][i], 0)
            r == 0 && continue          # agent i is gone by frame s
            base = late[r]
        end
        empty!(buf)
        inrange!(buf, tree, base, radius)
        zi = zt[i]
        for j in buf
            d = late[j] - base
            all(x -> abs(x) < W, d) || continue
            I = CartesianIndex(Tuple(round.(Int, d ./ dr) .+ (m + 1)))
            S[I, k] += zi * zs_[j]
            N[I, k] += 1
        end
    end
    return nothing
end
