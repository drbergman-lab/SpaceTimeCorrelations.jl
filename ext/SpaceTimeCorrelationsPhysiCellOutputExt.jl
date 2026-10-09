# PhysiCell output support, loaded when PhysiCellOutput is present.
#
# Adds a `spaceTimeCorrelation` method for a `PhysiCellSequence` (one output folder). It is a thin
# adapter: pick the cells, pull positions, values and IDs out of each snapshot's cell table, and
# call the array method. No second PhysiCell reader; PhysiCellOutput does the loading.

module SpaceTimeCorrelationsPhysiCellOutputExt

using SpaceTimeCorrelations
using PhysiCellOutput: PhysiCellSequence, loadCells!
using Tables: Tables

"""
    spaceTimeCorrelation(sequence::PhysiCellSequence; value, cell_types = nothing,
                         include_dead = false, id = :ID, kwargs...) -> SpaceTimeCorrelationResult

Estimate the space-time correlation of a cell variable across the snapshots of one PhysiCell
output folder. Requires `using PhysiCellOutput`.

- `value`: the cell variable, e.g. `:nicd` or any custom data label.
- `cell_types`: cell type names to keep (a name or a collection of names); `nothing` keeps all.
- `include_dead`: keep cells PhysiCell marks as dead.
- `id`: the ID column, which selects the tissue frame. Pass `id = nothing` for the lab frame.

Positions are `position_1`, `position_2` (and `position_3` with `dims = 3`). The remaining
keywords (`lags`, `dr`, `rmax`, `dims`, `mincount`, `ntasks`) are those of the array method;
`lags` count snapshots, so a lag of `k` is `k` times the save interval. Cell data is loaded into
`sequence` if it has not been already. Snapshots must be equally spaced in time.

# Examples
```julia
using PhysiCellOutput, SpaceTimeCorrelations
sequence = PhysiCellSequence("output"; include_cells = true)
r = spaceTimeCorrelation(sequence; value = :nicd, cell_types = "rgc",
                         lags = 0:6:60, dr = 10.0, rmax = 100.0)
```
"""
function SpaceTimeCorrelations.spaceTimeCorrelation(sequence::PhysiCellSequence;
                                                    value::Union{Symbol,AbstractString},
                                                    cell_types = nothing,
                                                    include_dead::Bool = false,
                                                    id::Union{Symbol,AbstractString,Nothing} = :ID,
                                                    dims::Integer = 2, kwargs...)
    snapshots = sequence.snapshots
    isempty(snapshots) && throw(ArgumentError("the sequence has no snapshots"))
    1 <= dims <= 3 || throw(ArgumentError("dims must be 1, 2 or 3 for PhysiCell output, got $dims"))
    _checkspacing([s.time for s in snapshots])
    keep = _celltypes(cell_types, sequence)
    any(s -> isempty(s.cells), snapshots) && loadCells!(sequence)

    value = Symbol(value)
    id = id === nothing ? nothing : Symbol(id)
    poscols = ntuple(k -> Symbol("position_", k), dims)
    needed = (value, poscols..., :dead, :cell_type_name, something(id, :dead))

    positions = Vector{Vector{NTuple{dims,Float64}}}(undef, length(snapshots))
    values = Vector{Vector{Float64}}(undef, length(snapshots))
    ids = id === nothing ? nothing : Vector{Vector{Int}}(undef, length(snapshots))
    for (t, snapshot) in enumerate(snapshots)
        cols = Tables.columns(snapshot.cells)
        names = Tables.columnnames(cols)
        nrows = Tables.rowcount(cols)
        if nrows == 0
            positions[t] = NTuple{dims,Float64}[]
            values[t] = Float64[]
            ids === nothing || (ids[t] = Int[])
            continue
        end
        for name in needed
            name in names || throw(ArgumentError(
                "snapshot $(snapshot.index) has no cell column $(repr(name)); see `sequence.labels` " *
                "for the cell variables PhysiCell saved"))
        end
        rows = findall(i -> (include_dead || !Tables.getcolumn(cols, :dead)[i]) &&
                            (keep === nothing || Tables.getcolumn(cols, :cell_type_name)[i] in keep),
                       1:nrows)
        coords = map(c -> Tables.getcolumn(cols, c), poscols)
        positions[t] = [map(c -> Float64(c[i]), coords) for i in rows]
        values[t] = Float64.(Tables.getcolumn(cols, value)[rows])
        ids === nothing || (ids[t] = Int.(Tables.getcolumn(cols, id)[rows]))
    end
    return spaceTimeCorrelation(positions, values; ids, dims, kwargs...)
end

_celltypes(::Nothing, sequence) = nothing
_celltypes(name::Union{AbstractString,Symbol}, sequence) = _celltypes([name], sequence)
function _celltypes(names, sequence)
    keep = Set(String.(names))
    known = Set(values(sequence.cell_type_to_name_dict))
    unknown = setdiff(keep, known)
    isempty(unknown) || throw(ArgumentError(
        "unknown cell type(s) $(join(sort!(collect(unknown)), ", ")); this output has " *
        join(sort!(collect(known)), ", ")))
    return keep
end

function _checkspacing(times)
    length(times) < 3 && return
    dt = times[2] - times[1]
    for k in 2:length(times)-1
        isapprox(times[k+1] - times[k], dt; rtol = 1e-6, atol = 1e-9) || throw(ArgumentError(
            "snapshots are not equally spaced in time (interval $dt, then " *
            "$(times[k+1] - times[k]) after t = $(times[k])); lags count snapshots, so they " *
            "would mean different times"))
    end
end

end
