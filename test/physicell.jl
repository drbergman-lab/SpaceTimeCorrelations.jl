# The PhysiCellOutput extension. The PhysiCellModelManager extension is not tested here: it only
# resolves a `Simulation` to a `PhysiCellSequence` and forwards, and PCMM would bring a large
# dependency tree into the test environment. The fixture is the initial state and
# first four snapshots (cell files only) of a 2D PhysiCell cortex run with a Notch–Delta model.

import PhysiCellOutput
using Tables: Tables

const FIXTURE = joinpath(@__DIR__, "fixtures", "physicell_output")

# The same selection done by hand, through the table method.
function handtables(sequence, value, celltype)
    map(sequence.snapshots) do s
        c = Tables.columns(s.cells)
        rows = findall(i -> !Tables.getcolumn(c, :dead)[i] &&
                            Tables.getcolumn(c, :cell_type_name)[i] == celltype,
                       1:Tables.rowcount(c))
        (ID = Tables.getcolumn(c, :ID)[rows], x = Tables.getcolumn(c, :position_1)[rows],
         y = Tables.getcolumn(c, :position_2)[rows], v = Tables.getcolumn(c, value)[rows])
    end
end

@testset "PhysiCellOutput extension" begin
    sequence = PhysiCellOutput.PhysiCellSequence(FIXTURE)   # cells not loaded yet
    @test length(sequence.snapshots) == 4
    kw = (lags = 0:2, dr = 10.0, rmax = 50.0)

    r = spaceTimeCorrelation(sequence; value = :notch, cell_types = "rgc", kw...)
    @test r.frame === :tissue
    @test all(s -> !isempty(s.cells), sequence.snapshots)   # loaded on demand
    ref = spaceTimeCorrelation(handtables(sequence, :notch, "rgc"); value = :v,
                               position = (:x, :y), id = :ID, kw...)
    @test isequal(r.C, ref.C)
    @test r.counts == ref.counts

    # names as strings, several types, the lab frame
    lab = spaceTimeCorrelation(sequence; value = "notch", cell_types = ["rgc"], id = nothing, kw...)
    @test lab.frame === :lab
    @test lab.counts[:, :, 1] == r.counts[:, :, 1]
    both = spaceTimeCorrelation(sequence; value = :notch, cell_types = ["rgc", "apical"], kw...)
    @test sum(both.counts) > sum(r.counts)
    all_ = spaceTimeCorrelation(sequence; value = :notch, kw...)
    @test sum(all_.counts) >= sum(both.counts)

    @test_throws ArgumentError spaceTimeCorrelation(sequence; value = :notch, cell_types = "rcg", kw...)
    @test_throws ArgumentError spaceTimeCorrelation(sequence; value = :not_a_label, kw...)
    @test_throws ArgumentError spaceTimeCorrelation(sequence; value = :notch, dims = 4, kw...)

    ext = Base.get_extension(SpaceTimeCorrelations, :SpaceTimeCorrelationsPhysiCellOutputExt)
    @test ext._checkspacing([0.0, 10.0, 20.0, 30.0]) === nothing
    @test_throws ArgumentError ext._checkspacing([0.0, 10.0, 20.0, 40.0])
end
