# PhysiCell output

With PhysiCellOutput loaded, `spaceTimeCorrelation` takes a whole output folder as a
`PhysiCellSequence`; with PhysiCellModelManager loaded, it takes a `Simulation`.

!!! tierbrief
    Name the cell variable with `value` and pick cell types with `cell_types`. The tissue frame
    is the default, since PhysiCell always records cell IDs.

!!! tierfull
    Restricting to the cell types that carry the signal usually matters. A variable that is
    constant in most cell types, as many custom variables are, makes an all-cell map mostly a
    map of where the varying cell type sits.

    Dead cells are left out unless `include_dead = true`. A frame in which every selected cell
    is dead has no z-score and is listed in `skipped`, like any other degenerate frame.

    `lags` count snapshots: with output saved every 10 minutes, `lags = 0:6:60` is 0 to 10 hours
    in steps of one hour. The snapshots must be equally spaced in time.

```julia
using PhysiCellOutput, SpaceTimeCorrelations

sequence = PhysiCellSequence("output"; include_cells = true)
r = spaceTimeCorrelation(sequence; value = :notch, cell_types = "rgc",
                         lags = 0:6:60, dr = 10.0, rmax = 100.0)
```

```julia
using PhysiCellModelManager, SpaceTimeCorrelations

r = spaceTimeCorrelation(Simulation(1); value = :notch, cell_types = ["rgc", "apical"],
                         lags = 0:6:60, dr = 10.0, rmax = 100.0)
```

!!! tierdev
    Both live in `ext/`. `SpaceTimeCorrelationsPhysiCellOutputExt` pulls `position_1`…,
    the value and `ID` out of each snapshot's cell table through Tables.jl (no DataFrames
    dependency) and calls the array method. `SpaceTimeCorrelationsPhysiCellModelManagerExt` only
    builds the `PhysiCellSequence` for a `Simulation` and forwards. The PhysiCellOutput extension
    is tested against `test/fixtures/physicell_output` (four snapshots of a real run, cell files
    only); the PCMM one is not tested here, because PCMM would bring a large dependency tree into
    the test environment.

!!! tierjournal "2026-10-08 — PhysiCell support as extensions, dead cells out by default"
    Decided: PhysiCellOutput and PhysiCellModelManager are weak dependencies, each reached through
    a thin extension over the array method, with a `cell_types` filter.
    Decided: exclude dead cells by default, matching PhysiCellOutput's `cellDataSequence`. In the
    first test run (a cortex model) every rgc cell was marked dead for 51 frames before removal;
    including them would correlate the values of cells that have stopped signalling.
    Rejected: testing the PCMM extension here. PCMM's dependencies conflict with the
    GeoStatsFunctions cross-check (Optim 1 against Optim 2), and the extension is one line.
