# Rayleigh–Ritz upper bounds, sparse ring blocks, partial spectra and coupled seed defects.

using CertifiedGaps.RitzBounds, CertifiedGaps.VerifiedSpectra, CertifiedGaps.RingModels
using CertifiedGaps.Seeds: coupled_seed_defects, seed_defects
using SparseArrays, Random
using CertifiedGaps.ExactArithmetic: RationalBig

@testset "Rayleigh–Ritz upper bounds" begin
    rng = MersenneTwister(7)
    n = 200
    A = sprandn(rng, n, n, 0.02)
    M = A + transpose(A) + spdiagm(0 => 4 .* randn(rng, n))
    @test issymmetric(M)
    # Reference: verified enclosures; every true eigenvalue is at least `exact` (lower ends).
    ref = verified_eigen(Matrix(M))
    exact = RationalBig.(ref.values) .- ref.radius
    F = eigen(Symmetric(Matrix(M)))
    for k in (1, 5, 20)
        r = ritz_upper_bounds(M, F.vectors[:, 1:k])
        @test length(r.upper) == k && issorted(r.upper)
        @test all(r.upper .>= exact[1:k]) && all(r.upper_float .>= r.upper)
        @test maximum(Float64.(r.upper - exact[1:k])) < 1e-9
        # A perturbed (non-orthonormal) basis still gives valid, slightly weaker bounds.
        V = F.vectors[:, 1:k] + 1e-6 * randn(rng, n, k)
        p = ritz_upper_bounds(M, V)
        @test all(p.upper .>= exact[1:k]) && p.orthogonality > 1e-7
        # A random orthonormalized basis gives valid but poor bounds.
        q = ritz_upper_bounds(M, Matrix(qr(randn(rng, n, k)).Q))
        @test all(q.upper .>= exact[1:k])
    end
    # Dense input and a matrix error: the bound moves up by exactly the error.
    r0 = ritz_upper_bounds(Matrix(M), F.vectors[:, 1:3])
    r1 = ritz_upper_bounds(Matrix(M), F.vectors[:, 1:3]; matrix_error = 1 // 1000)
    @test r1.upper == r0.upper .+ 1 // 1000
    # Shifting M changes nothing essential (the shift c follows Gershgorin).
    s = ritz_upper_bounds(M + 100I, F.vectors[:, 1:3])
    @test maximum(Float64.(s.upper .- 100 .- exact[1:3])) < 1e-9
    @test_throws ArgumentError ritz_upper_bounds(M + sparse([1], [2], [1e-3], n, n), F.vectors)
    @test_throws ErrorException ritz_upper_bounds(M, fill(1.0, n, 2) / sqrt(n))   # rank one
    # Lanczos on a large sparse tridiagonal Toeplitz matrix, exact spectrum cos(kπ/(m+1)).
    m = 1000
    T = spdiagm(1 => fill(-0.5, m - 1), -1 => fill(-0.5, m - 1))
    truth = sort([-cospi(big(k) / (m + 1)) for k in 1:m])
    b = lowest_upper_bounds(T, 6; tol = 1e-12, krylovdim = 120, maxiter = 1000,
        dense_max = 100)
    @test b.method == "ritz" && all(b.upper .>= truth[1:6])
    @test b.upper[1] - truth[1] < 1e-8
    # Small blocks are enclosed completely (also complex Hermitian ones).
    d = lowest_upper_bounds(Matrix(T[1:50, 1:50]) .+ 0im, 3)
    @test d.method == "dense" && length(d.upper) == 50
    @test all(d.upper .>= sort([-cospi(big(k) / 51) for k in 1:50]))
    @test_throws ArgumentError lowest_upper_bounds(sparse(T .+ 0im), 3; dense_max = 10)
end

@testset "sparse ring blocks and partial spectra" begin
    one = [1 0; 0 1]
    rotated(J) = J * (-kron(SX, SX) - kron(ISY, ISY) - kron(SZ, SZ))
    swap = [1 0 0 0; 0 0 1 0; 0 1 0 0; 0 0 0 1]
    dimer = RingModel(4, kron(one, rotated(1 // 2), one), rotated(3 // 2);
        charge = [0, 1, -1, 0], flip = kron(2SX, 2SX), reflection = swap * kron(2SX, 2SX))
    for L in (3, 4, 5)
        dense, sparse_blocks = Any[], Any[]
        ring_blocks((M, ε, mult, label) -> push!(dense, (M, ε, mult)), dimer, L)
        ring_blocks((M, ε, mult, label) -> push!(sparse_blocks, (M, ε, mult)), dimer, L;
            sparse = true)
        @test length(dense) == length(sparse_blocks)
        @test all(b[1] isa SparseMatrixCSC && Matrix(b[1]) == a[1] && a[2:3] == b[2:3]
                  for (a, b) in zip(dense, sparse_blocks))
        # A charge selection covers exactly the chosen sectors (q and −q for a flip pair).
        dims = Int[]
        ring_blocks((M, ε, mult, label) -> push!(dims, size(M, 1) * mult), dimer, L;
            exact = false, sparse = true, charges = [1])
        @test sum(dims) == 2 * sum(binomial(2L, L + 1))
    end
    # Partial spectrum of the 12-spin ring with forced Lanczos against the complete spectrum.
    L = 6
    full = verified_ring_spectrum(dimer, L)
    part = partial_ring_spectrum(dimer, L; howmany = 8, dense_max = 40)
    @test any(b.label["method"] == "ritz" for b in part)
    @test sum(b.label["dimension"] * b.multiplicity for b in part) == 4^L
    complete = partial_ring_spectrum(dimer, L; howmany = 4^L, dense_max = 10^4)
    for β in (1 // 2, 3, 8)
        z = log_partition_enclosure(full, β)
        @test inf(log_partition_lower(part, β)) <= sup(z)
        @test inf(log_partition_lower(complete, β)) <= sup(z)
        @test sup(z) - inf(log_partition_lower(complete, β)) < 1e-9 * max(1, β)
    end
    # Low-temperature weight is captured by 8 levels per block.
    @test sup(log_partition_enclosure(full, 8)) - inf(log_partition_lower(part, 8)) < 1e-6
    # Only some charges: still a lower bound.
    some = partial_ring_spectrum(dimer, L; howmany = 8, charges = [0], dense_max = 40)
    @test inf(log_partition_lower(some, 2)) < inf(log_partition_lower(part, 2))
    mktempdir() do dir
        path = write_partial_spectrum(joinpath(dir, "p.json.gz"), part)
        back = read_partial_spectrum(path)
        @test [b.upper for b in back] == [b.upper for b in part]
        @test partial_spectrum_summary(back)["levels"] == partial_spectrum_summary(part)["levels"]
    end
    # Coupled seed defects: εₛ as in `seed_defects`; εₜ via Z_{2n}(2β)/Z_n(β)⁴ ≥ the true
    # thermal purity times Pₛ², so it is larger than the direct one.
    # Even n = 2, so Pₛ(n, β) ≤ 1 by reflection positivity.
    n, β = 2, 2
    small, large = verified_ring_spectrum(dimer, n), verified_ring_spectrum(dimer, 2n)
    zn, z2n, z2n2 = (log_partition_enclosure(s, b) for (s, b) in
                     ((small, β), (large, β), (large, 2β)))
    direct = seed_defects(zn, z2n, z2n2)
    all4 = partial_ring_spectrum(dimer, 2n; howmany = 4^(2n), dense_max = 10^4)
    coupled = coupled_seed_defects(zn, log_partition_lower(all4, β),
        log_partition_lower(all4, 2β))
    @test coupled.eps_s >= direct.eps_s && coupled.eps_s - direct.eps_s < 1e-8
    @test coupled.eps_t >= direct.eps_t
    ps = exp(mid(z2n) - 2mid(zn))
    @test isapprox(1 - Float64(coupled.eps_t), (1 - Float64(direct.eps_t)) * ps^2; rtol = 1e-7)
    # Partial spectra only weaken the bounds.
    some4 = partial_ring_spectrum(dimer, 2n; howmany = 2, dense_max = 10)
    weak = coupled_seed_defects(zn, log_partition_lower(some4, β), log_partition_lower(some4, 2β))
    @test weak.eps_s > coupled.eps_s && weak.eps_t > coupled.eps_t
end

@testset "coupled seeds in gap certificates" begin
    # Staggered-field Heisenberg deep in the admissible region (see "end-to-end gap
    # certificate"); the doubled ring of 8 sites gets Ritz bounds in two charge sectors only.
    m = rotated_staggered(6)
    rp = rp_certificate(m; factors = [SX, ISY, SZ])
    small = verified_ring_spectrum(m, 4)
    large = partial_ring_spectrum(m, 8; howmany = 6, charges = [0, 1], dense_max = 20)
    seed = coupled_seed_bounds(m, 4, 1; small, large)
    complete = seed_bounds(m, 4, 1)
    @test seed.eps_s >= complete.eps_s && seed.eps_t >= complete.eps_t
    cert = certify_gap(m; n0 = 4, β0 = 1, rp, seed)
    @test cert["seed"]["route"] == "coupled" && haskey(cert["log_partition"],
        "logZ_2n0_2beta0_lower")
    stored = parse_exact_json(JSON3.write(cert))
    @test check_certificate(stored; spectra = false)["gap_lower_bound"] == cert["gap_lower_bound"]
    @test check_certificate(stored; spectra = Dict(4 => small, 8 => large))["defect_t_upper"] <=
          seed.eps_t
    # Its defects cannot be recomputed from the certificate alone, nor from weaker spectra.
    @test_throws ErrorException check_certificate(stored)
    weaker = partial_ring_spectrum(m, 8; howmany = 2, charges = [0], dense_max = 20)
    @test_throws ErrorException check_certificate(stored; spectra = Dict(4 => small, 8 => weaker))
    # An odd n0 needs a positive transfer operator (Pₛ ≤ 1).
    @test_throws ArgumentError coupled_seed_bounds(m, 3, 1; small, large)
end
