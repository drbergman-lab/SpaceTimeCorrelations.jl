"""
    SpaceTimeCorrelations

Estimate the joint space-time correlation `C(Δx, Δy, Δt)` of a continuous per-agent value on
moving, tracked agents.

The estimator is [`spaceTimeCorrelation`](@ref), which takes either plain per-frame arrays or
one Tables.jl table per frame, and returns a [`SpaceTimeCorrelationResult`](@ref). With
PhysiCellOutput or PhysiCellModelManager loaded, it also reads PhysiCell output directly.
"""
module SpaceTimeCorrelations

using NearestNeighbors: KDTree, inrange!
using StaticArrays: SVector
using Statistics: mean, std
using Tables: Tables

export spaceTimeCorrelation, SpaceTimeCorrelationResult

include("result.jl")
include("estimator.jl")
include("tables.jl")

end
