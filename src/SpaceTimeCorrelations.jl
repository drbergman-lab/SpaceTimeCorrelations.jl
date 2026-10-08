"""
    SpaceTimeCorrelations

Estimate the joint space-time correlation `C(Δx, Δy, Δt)` of a continuous per-agent value on
moving, tracked agents.

The estimator is [`spaceTimeCorrelation`](@ref), which takes either plain per-frame arrays or
one Tables.jl table per frame, and returns a [`SpaceTimeCorrelationResult`](@ref).
[`noiseFrames`](@ref) and [`bumpFrames`](@ref) generate toy data with known ground truth, and
[`frameTables`](@ref) turns it into per-frame tables.
"""
module SpaceTimeCorrelations

using NearestNeighbors: KDTree, inrange!
using Random: AbstractRNG, default_rng, randperm
using StaticArrays: SVector
using Statistics: mean, std
using Tables: Tables

export spaceTimeCorrelation, SpaceTimeCorrelationResult
export noiseFrames, bumpFrames, frameTables

include("result.jl")
include("estimator.jl")
include("tables.jl")
include("synthetic.jl")

end
