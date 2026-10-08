# An independent reference for the lag-0 slice. On z-scored values, the semivariogram of a pair
# is ½(zᵢ² + zⱼ²) - zᵢzⱼ, whose mean is close to 1 - C. GeoStatsFunctions computes the isotropic
# empirical variogram of one frame; averaging C over the grid bins at the same distances should
# match 1 - γ. GeoStatsFunctions is a test-only dependency: it covers lag 0 only, since it pairs
# values at the same locations within one point set.

import GeoStatsFunctions, GeoTables, Meshes
using Statistics: mean, std

@testset "lag 0 matches GeoStatsFunctions variogram" begin
    d = bumpFrames(; n = 3000, L = 40.0, nframes = 1, σ = 1.0, density = 0.1, rng = Xoshiro(2))
    v = d.values[1]
    z = (v .- mean(v)) ./ std(v; corrected = false)
    gt = GeoTables.georef((z = z,), [Meshes.Point(p[1], p[2]) for p in d.positions[1]])
    maxlag, nlags = 4.0, 8
    γ = only(GeoStatsFunctions.variogram(gt, :z; maxlag, nlags).ordinates)

    r = spaceTimeCorrelation(d.positions, d.values; lags = [0], dr = 0.25, rmax = maxlag + 0.5)
    edges = range(0, maxlag; length = nlags + 1)
    for k in 1:nlags
        s = 0.0; n = 0
        for ix in eachindex(r.centers), iy in eachindex(r.centers)
            edges[k] < hypot(r.centers[ix], r.centers[iy]) <= edges[k+1] || continue
            s += r.C[ix, iy, 1] * r.counts[ix, iy, 1]
            n += r.counts[ix, iy, 1]
        end
        @test s / n ≈ 1 - γ[k] atol = 0.05
    end
end
