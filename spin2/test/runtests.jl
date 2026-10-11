using Test, LinearAlgebra, ITensors, ITensorMPS
using Spin2Gap
include(joinpath(@__DIR__, "..", "scripts", "screen.jl"))
BLAS.set_num_threads(1)

@testset "necessary seed window" begin
    @test !beta_window(0.1, 1.0)["window_nonempty"]
    @test beta_window(0.001, 1.0)["window_nonempty"]
    @test isnothing(beta_window(0.0, 1.0)["beta_max_from_spatial_purity"])
    @test_throws ArgumentError beta_window(0.1, 0.0)
    @test_throws ArgumentError beta_window(0.1, 1.0; purity = 1.0)
end

@testset "spin-2 normalization and real transfer generators" begin
    s = spin_matrices()
    @test s.sx*s.sy-s.sy*s.sx ≈ im*s.sz
    @test s.sx^2+s.sy^2+s.sz^2 ≈ 6I
    gs = real_generators()
    @test all(g -> g ≈ -transpose(g), gs)
    @test all(g -> opnorm(g) ≈ 2, gs)
    @test gs[1]*gs[2]-gs[2]*gs[1] ≈ -gs[3]
    @test -sum(g*g for g in gs) ≈ 6I
    bond = kron(s.sz, s.sz)+(kron(s.sp, s.sm)+kron(s.sm, s.sp))/2
    @test eigvals(Symmetric(bond)) ≈
          vcat(fill(-6.0, 1), fill(-5.0, 3), fill(-3.0, 5), fill(0.0, 7), fill(4.0, 9)) atol=1e-13
    @test eigvals(Symmetric(-sum(kron(g, g) for g in gs))) ≈ eigvals(Symmetric(bond)) atol=1e-13
end

@testset "MPO and independent Hilbert-space reference" begin
    for periodic in (false, true)
        n = 4
        sites = siteinds("S=2", n)
        H = hamiltonian(sites; periodic)
        dense = reshape(Array(prod(H), prime.(sites)..., dag.(sites)...), 5^n, 5^n)
        reference = exact_spectrum(n; periodic)
        @test eigvals(Symmetric(dense)) ≈ reference atol=2e-12
        @test length(reference) == 5^n
        @test log_partition(reference, 0) ≈ n*log(5)
        if periodic
            # H4=(S1+S3).(S2+S4), hence E0=-20 and first gap=1.
            @test reference[1] ≈ -20
            @test reference[2]-reference[1] ≈ 1
        end
    end
    @test_throws ArgumentError bonds(2)
    @test_throws ArgumentError bonds(5)
    @test_throws ArgumentError sector_hamiltonian(8, 0)
end

@testset "thermal normalization, periodic closure, and Trotter convergence" begin
    es = exact_spectrum(4)
    beta = 0.1
    a = thermal_curve(4, [0.0, beta]; dbeta = 0.05, maxdim = 625, cutoff = 1e-14, progress = false)
    b = thermal_curve(4, [beta]; dbeta = 0.025, maxdim = 625, cutoff = 1e-14, progress = false)
    truth = log_partition(es, beta)
    @test a[0.0]["logZ"] ≈ 4log(5)
    err1 = abs(a[beta]["logZ"]-truth)
    err2 = abs(b[beta]["logZ"]-truth)
    @test err2 < err1/3
    @test err2 < 1e-4
    # A single open bond has commuting evolution and a known full spectrum.
    e2 = exact_spectrum(2; periodic = false)
    c = thermal_curve(
        2,
        [0.2, 0.4];
        dbeta = 0.1,
        maxdim = 25,
        cutoff = 1e-14,
        periodic = false,
        progress = false,
    )
    @test c[0.2]["logZ"] ≈ log_partition(e2, 0.2) atol=1e-11
    @test c[0.4]["logZ"] ≈ log_partition(e2, 0.4) atol=1e-11
    @test_throws ArgumentError thermal_curve(4, [0.13]; dbeta = 0.1)
end

@testset "QN DMRG and excited-state penalty against full diagonalization" begin
    row = dmrg_levels(4; maxdim = 32, nsweeps = 6, cutoff = 1e-12, outputlevel = 0)
    @test row["levels"]["ground_sz0"]["energy"] ≈ -20 atol=1e-8
    @test row["sz1_gap_estimate"] ≈ 1 atol=1e-8
    @test row["sz0_gap_estimate"] ≈ 1 atol=1e-8
    @test row["levels"]["excited_sz0"]["ground_overlap"] < 1e-6
    @test all(abs(v["variance"]) < 1e-7 for v in values(row["levels"]))
    @test row["certified"] == false
end
