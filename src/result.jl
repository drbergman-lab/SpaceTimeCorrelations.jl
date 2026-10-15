"""
    SpaceTimeCorrelationResult

The output of [`spaceTimeCorrelation`](@ref).

# Fields
- `C::Array{Float64}`: the correlation, indexed `C[ix, iy, k]` in 2D (`C[ix, iy, iz, k]` in 3D).
  Spatial indices run over `centers`, `k` over `lags`. Bins with fewer than `min_pairs` pairs are
  `NaN`.
- `counts::Array{Int}`: the number of pairs in each bin, same shape as `C`.
- `centers::Vector{Float64}`: the bin centers along every spatial axis. The grid is square, so
  one vector serves all axes; `centers[n + 1] == 0` where `length(centers) == 2n + 1`.
- `lags::Vector{Int}`: the time lags, in frames, in the order given.
- `dr::Float64`: the bin width.
- `frame::Symbol`: `:lab` or `:tissue`, the reference frame displacements are measured in.
- `skipped::Vector{Int}`: frames left out because their z-score is undefined (fewer than two
  agents, or zero variance).

# Examples
```jldoctest
julia> positions = [[(0.0, 0.0), (1.0, 0.0), (0.0, 1.0)], [(0.0, 0.0), (1.0, 0.0), (0.0, 1.0)]];

julia> values = [[1.0, 2.0, 3.0], [3.0, 2.0, 1.0]];

julia> r = spaceTimeCorrelation(positions, values; lags = 0:1, dr = 0.5, rmax = 2);

julia> size(r.C), r.centers[5], r.frame
((9, 9, 2), 0.0, :lab)
```
"""
struct SpaceTimeCorrelationResult{N}
    C::Array{Float64,N}
    counts::Array{Int,N}
    centers::Vector{Float64}
    lags::Vector{Int}
    dr::Float64
    frame::Symbol
    skipped::Vector{Int}
end

function Base.show(io::IO, ::MIME"text/plain", r::SpaceTimeCorrelationResult{N}) where {N}
    nb = length(r.centers)
    println(io, "SpaceTimeCorrelationResult ($(r.frame) frame, $(N - 1)D)")
    println(io, "  bins:  ", join(fill(nb, N - 1), " × "), ", width ", r.dr,
            ", centers ", first(r.centers), " … ", last(r.centers))
    println(io, "  lags:  ", r.lags)
    print(io, "  pairs: ", sum(r.counts))
    isempty(r.skipped) || print(io, "\n  skipped frames: ", r.skipped)
end
