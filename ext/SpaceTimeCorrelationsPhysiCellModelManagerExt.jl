# PhysiCellModelManager support, loaded when PhysiCellModelManager and PhysiCellOutput are present.
#
# A thin adapter: PCMM resolves a `Simulation` to its output folder and loads it as a
# `PhysiCellSequence`; the PhysiCellOutput extension does the rest.

module SpaceTimeCorrelationsPhysiCellModelManagerExt

using SpaceTimeCorrelations
using PhysiCellModelManager: Simulation
using PhysiCellOutput: PhysiCellSequence

"""
    spaceTimeCorrelation(simulation::Simulation; value, cell_types = nothing, kwargs...)
        -> SpaceTimeCorrelationResult

Estimate the space-time correlation of a cell variable across one PhysiCellModelManager
simulation's output. Requires `using PhysiCellModelManager`.

Loads the simulation's output as a `PhysiCellSequence` with cell data and forwards every keyword
to [`spaceTimeCorrelation(::PhysiCellSequence)`](@ref
SpaceTimeCorrelations.spaceTimeCorrelation(::PhysiCellOutput.PhysiCellSequence)).

# Examples
```julia
using PhysiCellModelManager, SpaceTimeCorrelations
r = spaceTimeCorrelation(Simulation(1); value = :nicd, cell_types = "rgc",
                         lags = 0:6:60, dr = 10.0, rmax = 100.0)
```
"""
SpaceTimeCorrelations.spaceTimeCorrelation(simulation::Simulation; kwargs...) =
    spaceTimeCorrelation(PhysiCellSequence(simulation; include_cells = true); kwargs...)

end
