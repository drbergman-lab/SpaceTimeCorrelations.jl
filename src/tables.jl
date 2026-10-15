"""
    spaceTimeCorrelation(frames; value, position = :position, id = nothing, kwargs...)
        -> SpaceTimeCorrelationResult

Estimate the space-time correlation from a vector of tables, one per frame, in time order.

Each element of `frames` is any Tables.jl table (a `DataFrame`, a `NamedTuple` of vectors, a
CSV file, …) with one row per agent. Column names are passed as `Symbol`s:

- `value`: the column holding the per-agent value. Required.
- `position`: either one column whose entries are points `(x, y[, z])`, or a tuple of
  coordinate columns such as `(:x, :y)`.
- `id`: the column of persistent agent IDs. Passing it selects the tissue frame.

The remaining keywords (`lags`, `dr`, `rmax`, `dims`, `mincount`, `ntasks`) are those of the
array method, [`spaceTimeCorrelation(positions, values)`](@ref
spaceTimeCorrelation(::AbstractVector, ::AbstractVector)).

# Examples
```jldoctest
julia> frames = [(id = [1, 2, 3], x = [0.0, 1.0, 0.0], y = [0.0, 0.0, 1.0], value = [1.0, 2.0, 3.0]),
                 (id = [1, 2, 3], x = [0.0, 1.0, 0.0], y = [0.0, 0.0, 1.0], value = [3.0, 2.0, 1.0])];

julia> r = spaceTimeCorrelation(frames; value = :value, position = (:x, :y), id = :id,
                                lags = 0:1, dr = 1, rmax = 2);

julia> r.frame
:tissue
```
"""
function spaceTimeCorrelation(frames::AbstractVector; value::Symbol,
                              position::Union{Symbol,Tuple{Vararg{Symbol}}} = :position,
                              id::Union{Symbol,Nothing} = nothing, kwargs...)
    cols = map(Tables.columns, frames)
    positions = [_positioncolumn(c, position, t) for (t, c) in enumerate(cols)]
    values = [_column(c, value, t) for (t, c) in enumerate(cols)]
    ids = id === nothing ? nothing : [_column(c, id, t) for (t, c) in enumerate(cols)]
    return spaceTimeCorrelation(positions, values; ids, kwargs...)
end

function _column(cols, name::Symbol, t)
    name in Tables.columnnames(cols) || throw(ArgumentError(
        "frame $t has no column $(repr(name)); its columns are $(Tuple(Tables.columnnames(cols)))"))
    return Tables.getcolumn(cols, name)
end

_positioncolumn(cols, name::Symbol, t) = _column(cols, name, t)

function _positioncolumn(cols, names::Tuple{Vararg{Symbol}}, t)
    isempty(names) && throw(ArgumentError("position must name at least one column"))
    coords = map(n -> _column(cols, n, t), names)
    return [map(c -> c[i], coords) for i in eachindex(first(coords))]
end
