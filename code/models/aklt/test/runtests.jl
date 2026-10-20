using Test, LinearAlgebra, ITensors, ITensorMPS
using AKLTPilot
using CertifiedGaps.Spins: spin_matrices
BLAS.set_num_threads(1)

@testset "honeycomb quotient and coloring" begin
    for (nx, ny) in ((2, 2), (4, 2), (4, 4))
        graph = honeycomb(nx, ny)
        @test graph.nsites==2nx*ny
        @test length(graph.bonds)==3nx*ny
        for group in graph.groups
            @test sort(vcat([[a, b] for (a, b) in group]...))==collect(1:graph.nsites)
        end
        @test all(graph.coordinates[a][3]!=graph.coordinates[b][3] for (a, b) in graph.bonds)
    end
    @test_throws ArgumentError honeycomb(1, 2)
end

@testset "projector, normalization and reflection witness" begin
    p = bond_projector()
    s = spin_matrices(3//2)
    x = kron(s.sz, s.sz)+(kron(s.sp, s.sm)+kron(s.sm, s.sp))/2
    @test p ≈ (x+15I/4)*(x+11I/4)*(x+3I/4)/90 atol=1e-14
    @test p*p ≈ p atol=1e-14
    @test tr(p) ≈ 7
    @test s.sx^2+s.sy^2+s.sz^2 ≈ 15I/4
    q = s.sz^2-5I/4
    @test tr(p*kron(q, q)) ≈ 7/5
    for t in (0.1, 1.0)
        @test exp(-t*Symmetric(p)) ≈ I+expm1(-t)*p atol=1e-14
    end
end

@testset "thermal full trace and nonlocal gates" begin
    # A single bond has nine zero eigenvalues and seven eigenvalues equal to one.
    curves = thermal_curve(
        2,
        [(1, 2)],
        [0.0, 0.25, 0.5];
        dbeta = 0.125,
        maxdim = 16,
        cutoff = 1e-14,
        progress = false,
    )
    for beta in (0.0, 0.25, 0.5)
        @test curves[beta]["logZ"] ≈ log(9+7exp(-beta)) atol=2e-12
    end
    # Four-site ring uses the nonlocal closing bond and noncommuting projectors.
    bonds = [(1, 2), (2, 3), (3, 4), (1, 4)]
    exact = full_spectrum(4, bonds; progress = false).energies
    @test length(exact)==256
    @test log_partition(exact, 0) ≈ 4log(4)
    beta = 0.25
    a = thermal_curve(
        4,
        bonds,
        [beta];
        dbeta = 0.125,
        maxdim = 256,
        cutoff = 1e-14,
        progress = false,
    )
    b = thermal_curve(
        4,
        bonds,
        [beta];
        dbeta = 0.0625,
        maxdim = 256,
        cutoff = 1e-14,
        progress = false,
    )
    truth = log_partition(exact, beta)
    ea, eb = abs(a[beta]["logZ"]-truth), abs(b[beta]["logZ"]-truth)
    @test ea>1e-10
    @test eb<ea/3
    @test eb<1e-5
    # Keeping moved sites between gates changes the compression path, not exact evolution.
    c = thermal_curve(
        4,
        bonds,
        [beta];
        dbeta = 0.0625,
        maxdim = 256,
        cutoff = 1e-14,
        progress = false,
        move_sites_back_between_gates = false,
    )
    # Different SVD sequences accumulate roundoff and cutoff-level errors.
    @test c[beta]["logZ"] ≈ b[beta]["logZ"] atol=1e-10
    @test_throws ArgumentError thermal_curve(2, [(1, 2)], [0.2]; dbeta = 0.125)
end

@testset "diagnostic does not conceal invalid normalization" begin
    good = Dict(1.0=>Dict("logZ"=>0.01), 2.0=>Dict("logZ"=>0.001))
    bad = Dict(1.0=>Dict("logZ"=>-0.01), 2.0=>Dict("logZ"=>-0.001))
    @test isnothing(diagnostics(bad, good, good, 1.0)["C"])
    @test diagnostics(good, good, good, 1.0)["G_x"]<0
end
