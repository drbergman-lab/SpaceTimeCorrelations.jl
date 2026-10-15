using Test
using SpaceTimeCorrelations
using Random
using StaticArrays

include("synthetic.jl")   # toy-data generators with known answers

# Location of the correlation peak at lag slot k: the C-weighted centroid of the bins within half
# of the maximum, among bins with enough pairs to trust. More stable than argmax on a broad peak.
function peak(r, k; mincount = 200)
    C = r.C[:, :, k]
    ok = (r.counts[:, :, k] .>= mincount) .& isfinite.(C)
    top = maximum(C[ok])
    w = 0.0; x = 0.0; y = 0.0
    for I in CartesianIndices(C)
        ok[I] && C[I] >= top / 2 || continue
        w += C[I]; x += C[I] * r.centers[I[1]]; y += C[I] * r.centers[I[2]]
    end
    return (x / w, y / w)
end

zerobin(r) = (length(r.centers) + 1) ÷ 2
binat(r, v) = findfirst(c -> isapprox(c, v; atol = 1e-9), r.centers)
near(a, b; atol = 0.3) = all(isapprox.(a, b; atol))

const VELOCITY_KW = (n = 4000, L = 40.0, nframes = 6, σ = 1.0, density = 0.1)
const GRID = (lags = 0:3, dr = 0.5, rmax = 5)

@testset "SpaceTimeCorrelations" begin

@testset "grid" begin
    positions = [[(0.0, 0.0), (5.0, 0.0)], [(1.0, 0.0), (6.0, 0.0)]]
    values = [[1.0, -1.0], [1.0, -1.0]]
    r = spaceTimeCorrelation(positions, values; lags = [1], dr = 1, rmax = 2)
    @test r.centers == [-2.0, -1.0, 0.0, 1.0, 2.0]
    @test r.centers[zerobin(r)] === 0.0
    @test size(r.C) == (5, 5, 1)

    # rmax rounds out to whole bins
    r = spaceTimeCorrelation(positions, values; lags = [0], dr = 1, rmax = 1.2)
    @test r.centers == -2.0:2.0

    # a displacement just inside half a bin lands in the zero bin, one just past it does not
    p = [[(0.0, 0.0), (10.0, 0.0)], [(0.49, 0.0), (10.51, 0.0)]]
    r = spaceTimeCorrelation(p, values; lags = [1], dr = 1, rmax = 2)
    @test r.counts[3, 3, 1] == 1
    @test r.counts[4, 3, 1] == 1

    # corner bins are searched: the grid is the full square, not the inscribed disk
    p = [[(0.0, 0.0), (10.0, 0.0)], [(2.0, 2.0), (12.0, -2.0)]]
    r = spaceTimeCorrelation(p, values; lags = [1], dr = 1, rmax = 2)
    @test r.counts[5, 5, 1] == 1
    @test r.counts[5, 1, 1] == 1
end

@testset "sign convention" begin
    # The pattern (agent with value +1) moves by +1 in x: the peak is at positive Δx.
    positions = [[(0.0, 0.0), (5.0, 0.0)], [(1.0, 0.0), (6.0, 0.0)]]
    values = [[1.0, -1.0], [1.0, -1.0]]
    r = spaceTimeCorrelation(positions, values; lags = [1], dr = 1, rmax = 2)
    @test r.C[binat(r, 1.0), zerobin(r), 1] == 1.0
    @test r.counts[binat(r, 1.0), zerobin(r), 1] == 2
    @test isnan(r.C[binat(r, -1.0), zerobin(r), 1])
    @test r.frame === :lab

    # Swap the frames: the same motion run backwards puts it at -Δx.
    r = spaceTimeCorrelation(reverse(positions), reverse(values); lags = [1], dr = 1, rmax = 2)
    @test r.C[binat(r, -1.0), zerobin(r), 1] == 1.0
end

@testset "symmetry at lag 0" begin
    d = bumpFrames(; n = 1500, L = 30.0, nframes = 2, rng = Xoshiro(11))
    for ids in (nothing, d.ids)
        r = spaceTimeCorrelation(d.positions, d.values; ids, lags = [0, 1], dr = 0.5, rmax = 4)
        C = r.C[:, :, 1]
        @test isequal(r.counts[:, :, 1], reverse(r.counts[:, :, 1]))
        @test all(isapprox.(C, reverse(C); atol = 1e-12) .| (isnan.(C) .& isnan.(reverse(C))))
    end
end

@testset "1. noise baseline" begin
    d = noiseFrames(; n = 400, nframes = 5, rng = Xoshiro(1))
    r = spaceTimeCorrelation(d.positions, d.values; lags = 0:2, dr = 0.1, rmax = 3)
    m = zerobin(r)
    @test r.C[m, m, 1] ≈ 1 atol = 0.05
    for k in 1:3
        C = r.C[:, :, k]
        N = r.counts[:, :, k]
        mask = trues(size(C)); k == 1 && (mask[m, m] = false)
        mask .&= N .> 0
        @test abs(sum(C[mask] .* N[mask]) / sum(N[mask])) < 0.02
    end
end

@testset "2. cluster size" begin
    σ = 1.5
    d = bumpFrames(; n = 4000, L = 40.0, nframes = 1, σ, density = 0.05, rng = Xoshiro(2))
    r = spaceTimeCorrelation(d.positions, d.values; lags = [0], dr = 0.5, rmax = 5)
    m = zerobin(r)
    for x in 0.5:0.5:4.0
        i = binat(r, x)
        # average the four axis directions, which share the same distance
        c = (r.C[m+(i-m), m, 1] + r.C[m-(i-m), m, 1] + r.C[m, m+(i-m), 1] + r.C[m, m-(i-m), 1]) / 4
        @test c ≈ exp(-x^2 / (4σ^2)) atol = 0.12
    end
end

@testset "3. translating pattern" begin
    v = (1.0, 0.0)
    d = bumpFrames(; VELOCITY_KW..., relative_velocity = v, rng = Xoshiro(3))
    lab = spaceTimeCorrelation(d.positions, d.values; GRID...)
    tissue = spaceTimeCorrelation(d.positions, d.values; ids = d.ids, GRID...)
    for (k, lag) in enumerate(GRID.lags)
        @test near(peak(lab, k), v .* lag)
    end
    @test isequal(lab.C, tissue.C)  # fixed agents: the two frames coincide
    @test tissue.frame === :tissue
end

@testset "4. pure advection" begin
    u = (1.0, 0.5)
    d = bumpFrames(; VELOCITY_KW..., flow = u, rng = Xoshiro(4))
    lab = spaceTimeCorrelation(d.positions, d.values; GRID...)
    tissue = spaceTimeCorrelation(d.positions, d.values; ids = d.ids, GRID...)
    for (k, lag) in enumerate(GRID.lags)
        @test near(peak(lab, k), u .* lag)
        @test near(peak(tissue, k), (0.0, 0.0))
    end
end

@testset "5. relative motion" begin
    u = (1.0, 0.0)
    v = (0.0, 1.0)
    d = bumpFrames(; VELOCITY_KW..., flow = u, relative_velocity = v, rng = Xoshiro(5))
    lab = spaceTimeCorrelation(d.positions, d.values; GRID...)
    tissue = spaceTimeCorrelation(d.positions, d.values; ids = d.ids, GRID...)
    for (k, lag) in enumerate(GRID.lags)
        @test near(peak(tissue, k), v .* lag)
        @test near(peak(lab, k), (u .+ v) .* lag)
    end
end

@testset "6. diagonal motion" begin
    v = (1.0, -0.5)
    d = bumpFrames(; VELOCITY_KW..., relative_velocity = v, rng = Xoshiro(6))
    r = spaceTimeCorrelation(d.positions, d.values; ids = d.ids, GRID...)
    for (k, lag) in enumerate(GRID.lags)
        @test near(peak(r, k), v .* lag)
    end
    # the peak bin itself is off both axes
    k = 3
    C = copy(r.C[:, :, k]); C[r.counts[:, :, k] .< 200] .= -Inf
    I = argmax(C)
    @test (r.centers[I[1]], r.centers[I[2]]) == (2.0, -1.0)
end

@testset "7. opposing flows" begin
    L = 40.0
    flow = x -> x[2] < L / 2 ? (1.0, 0.0) : (-1.0, 0.0)
    d = bumpFrames(; VELOCITY_KW..., L, flow, rng = Xoshiro(7))
    lab = spaceTimeCorrelation(d.positions, d.values; GRID...)
    tissue = spaceTimeCorrelation(d.positions, d.values; ids = d.ids, GRID...)
    m = zerobin(lab)
    for (k, lag) in enumerate(GRID.lags)
        @test near(peak(tissue, k), (0.0, 0.0))
        lag < 2 && continue
        # two offset peaks, each higher than the centre between them
        plus = lab.C[binat(lab, lag), m, k]
        minus = lab.C[binat(lab, -lag), m, k]
        @test plus > 0.5 && minus > 0.5
        @test lab.C[m, m, k] < min(plus, minus) - 0.2
    end
end

@testset "8. turnover" begin
    base = bumpFrames(; VELOCITY_KW..., rng = Xoshiro(8))
    churn = bumpFrames(; VELOCITY_KW..., turnover = 0.2, rng = Xoshiro(8))
    @test churn.ids[2] != churn.ids[1]
    lags = 0:5
    r0 = spaceTimeCorrelation(base.positions, base.values; ids = base.ids, lags, dr = 0.5, rmax = 3)
    r = spaceTimeCorrelation(churn.positions, churn.values; ids = churn.ids, lags, dr = 0.5, rmax = 3)
    pairs = vec(sum(r.counts; dims = (1, 2)))
    pairs0 = vec(sum(r0.counts; dims = (1, 2)))
    @test pairs[1] == pairs0[1]
    @test all(pairs[2:end] .< pairs0[2:end])
    # counts per frame pair fall with lag as fewer base agents survive
    perpair = pairs ./ (6 .- collect(lags))
    @test issorted(perpair; rev = true)
    @test near(peak(r, 4), (0.0, 0.0))
end

@testset "degenerate frames are skipped" begin
    d = noiseFrames(; n = 50, nframes = 4, rng = Xoshiro(9))
    values = copy(d.values)
    values[2] = fill(3.0, 50)                  # zero variance
    positions = copy(d.positions)
    positions[4] = positions[4][1:1]           # one agent
    values[4] = values[4][1:1]
    r = spaceTimeCorrelation(positions, values; lags = 0:1, dr = 1, rmax = 2)
    @test r.skipped == [2, 4]
    # only frames 1 and 3 remain: two frame pairs at lag 0, none at lag 1
    @test sum(r.counts[:, :, 2]) == 0
    @test all(isnan, r.C[:, :, 2])
    @test sum(r.counts[:, :, 1]) > 0
    @test occursin("skipped frames: [2, 4]", sprint(show, MIME"text/plain"(), r))
end

@testset "mincount" begin
    d = noiseFrames(; n = 200, nframes = 3, rng = Xoshiro(10))
    r1 = spaceTimeCorrelation(d.positions, d.values; lags = 0:1, dr = 0.5, rmax = 2)
    r = spaceTimeCorrelation(d.positions, d.values; lags = 0:1, dr = 0.5, rmax = 2, mincount = 20)
    @test r.counts == r1.counts
    @test all(isnan.(r.C) .== ((r.counts .< 20) .| isnan.(r1.C)))
end

@testset "ntasks" begin
    d = bumpFrames(; n = 1000, L = 25.0, nframes = 5, flow = (0.5, 0.0), rng = Xoshiro(12))
    r1 = spaceTimeCorrelation(d.positions, d.values; ids = d.ids, lags = 0:3, dr = 0.5, rmax = 3,
                              ntasks = 1)
    r3 = spaceTimeCorrelation(d.positions, d.values; ids = d.ids, lags = 0:3, dr = 0.5, rmax = 3,
                              ntasks = 3)
    @test r1.counts == r3.counts
    @test all(isapprox.(r1.C, r3.C; atol = 1e-12) .| (isnan.(r1.C) .& isnan.(r3.C)))
end

@testset "input forms" begin
    d = bumpFrames(; n = 500, L = 20.0, nframes = 3, flow = (0.5, 0.0), rng = Xoshiro(13))
    kw = (lags = 0:2, dr = 0.5, rmax = 2)
    ref = spaceTimeCorrelation(d.positions, d.values; ids = d.ids, kw...)

    # matrix positions, one column per agent
    mats = [reduce(hcat, Vector.(p)) for p in d.positions]
    @test isequal(spaceTimeCorrelation(mats, d.values; ids = d.ids, kw...).C, ref.C)

    # 3-vectors with a z coordinate: dims = 2 uses x and y only
    p3 = [[(q[1], q[2], 7.0) for q in p] for p in d.positions]
    @test isequal(spaceTimeCorrelation(p3, d.values; ids = d.ids, kw...).C, ref.C)

    # table method, with separate coordinate columns and with one column of points
    tables = frameTables(d)
    rt = spaceTimeCorrelation(tables; value = :value, position = (:x, :y), id = :id, kw...)
    @test isequal(rt.C, ref.C)
    pointtables = [(cell = t.id, pos = collect(zip(t.x, t.y)), NICD = t.value) for t in tables]
    rp = spaceTimeCorrelation(pointtables; value = :NICD, position = :pos, id = :cell, kw...)
    @test isequal(rp.C, ref.C)
    @test spaceTimeCorrelation(tables; value = :value, position = (:x, :y), kw...).frame === :lab

    # 3D falls out of the same code
    r3 = spaceTimeCorrelation(p3, d.values; ids = d.ids, dims = 3, kw...)
    @test size(r3.C) == (9, 9, 9, 3)
    @test sum(r3.counts) == sum(ref.counts)
end

@testset "argument errors" begin
    d = noiseFrames(; n = 20, nframes = 2, rng = Xoshiro(14))
    kw = (lags = 0:1, dr = 1, rmax = 2)
    @test_throws ArgumentError spaceTimeCorrelation(d.positions, d.values[1:1]; kw...)
    @test_throws ArgumentError spaceTimeCorrelation(d.positions, d.values; ids = d.ids[1:1], kw...)
    @test_throws ArgumentError spaceTimeCorrelation(d.positions, d.values; lags = [-1], dr = 1, rmax = 2)
    @test_throws ArgumentError spaceTimeCorrelation(d.positions, d.values; lags = [1, 1], dr = 1, rmax = 2)
    @test_throws ArgumentError spaceTimeCorrelation(d.positions, d.values; lags = 0:1, dr = 0, rmax = 2)
    @test_throws ArgumentError spaceTimeCorrelation(d.positions, d.values; lags = 0:1, dr = 1, rmax = 0)
    @test_throws ArgumentError spaceTimeCorrelation(d.positions, d.values; dims = 3, kw...)
    bad = copy(d.values); bad[1] = [NaN; bad[1][2:end]]
    @test_throws ArgumentError spaceTimeCorrelation(d.positions, bad; kw...)
    dup = [copy(i) for i in d.ids]; dup[2][2] = dup[2][1]
    @test_throws ArgumentError spaceTimeCorrelation(d.positions, d.values; ids = dup, kw...)
    @test_throws ArgumentError spaceTimeCorrelation(frameTables(d); value = :NICD, kw...)
    @test_throws ArgumentError bumpFrames(; turnover = 2)
end

include("crosscheck_geostats.jl")
include("physicell.jl")

end
