# Fast checks of the seed-scale model functions and of the exact admissibility bisection.
using Test
include(joinpath(@__DIR__, "..", "experiments", "seed_scale.jl"))
using .SeedScale

@testset "magnon-gas model" begin
    @test SeedScale.dispersion(π) ≈ SeedScale.GAP
    @test SeedScale.dispersion(0.0) ≈ sqrt(SeedScale.GAP^2 + (SeedScale.VELOCITY * π)^2)
    @test magnon_sum(256, 64.0) < magnon_sum(256, 32.0)          # decreasing in β
    @test magnon_sum(512, 64.0) > magnon_sum(256, 64.0)          # more modes on the larger ring
    ls, lt = log_purities(256, 64.0, 4.7e-4)
    @test ls <= 0 && lt <= 0
    εs, εt = SeedScale.defects(256, 64.0, 4.7e-4)
    @test 0 < εs < 0.1 && 0 < εt < 0.05
end

@testset "ring closure from the certified trial energies" begin
    e = trial_energies()
    @test all(haskey(e, N) for N in (128, 256, 512, 1024))
    @test ring_closure(e, 128) ≈ 0.010644 atol = 1e-6
    @test ring_closure(e, 256) ≈ 4.694e-4 atol = 1e-6
    m = delta_model(e)
    @test m.δ(128) ≈ ring_closure(e, 128)
    @test m.δ(256) ≈ ring_closure(e, 256)
    @test 35 < m.xi < 50
    @test delta_model(e; xi = 49.49).δ(256) > m.δ(256)
end

@testset "exact admissible region" begin
    @test SeedScale.admissible(1 // 10, 1 // 10)
    @test !SeedScale.admissible(1 // 5, 1 // 5)
    r = max_defect(1 // 10^8; vary = :spatial, tol = 1 // 10^3)
    @test 1 // 5 < r.accepted < r.rejected < 1 // 4
    @test SeedScale.admissible(r.accepted, 1 // 10^8) &&
          !SeedScale.admissible(r.rejected, 1 // 10^8)
    @test SeedScale.squaring_budget(28, 0.02) ≈ 0.02 / 2^28
end

@testset "windows" begin
    w = seed_window(256, 40.0:8.0:120.0, 4.7e-4)
    @test w.βmin !== nothing && w.best_bound > 0.02
    @test seed_window(128, 40.0:8.0:120.0, 0.0106).βmin === nothing
end
