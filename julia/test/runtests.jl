using CertifiedGaps, IntervalArithmetic, LinearAlgebra, Test
using CertifiedGaps.ExactArithmetic, CertifiedGaps.Results, CertifiedGaps.Spins
using CertifiedGaps.ReflectionPositivity, CertifiedGaps.ExactDiagonalization
using CertifiedGaps.Models, CertifiedGaps.Seeds
using ITensors, ITensorMPS, CertifiedGaps.Purification

@testset "CertifiedGaps" begin
    @testset "recurrence constants (survey §2.6)" begin
        @test fixed_point_bracket_holds()
        εstar = diagonal_fixed_point()
        @test 0.12731 < εstar < 0.12732
        @test diagonal_defect(εstar) ≈ εstar rtol = 1e-30
        # On the diagonal below ε⋆ the map contracts, above it it does not (Theorem 2.17).
        @test all(doubling_map(0.12, 0.12) .< 0.12)
        @test all(doubling_map(0.13, 0.13) .>= 0.13)
        Λ, ok = setprecision(BigFloat, 256) do
            majorant_loss(interval(BigFloat, 3) / 97, interval(BigFloat, 1) / 3; J = 40)
        end
        @test ok && 0.129021 < sup(Λ) < 0.129022
    end

    @testset "certified bound" begin
        c = certified_gap_bound(1 // 100, 1 // 100, 1)
        @test c.j_admissible == 0
        @test c.bound <= gap_bound(1 // 100, 1 // 100, 1) <= c.upper
        @test 3.3642 < c.bound < 3.3643                      # survey: 3.3643 (notes: log 20)
        @test certified_gap_bound(1 // 100, 1 // 100, 784).bound > log(20) / 784
        @test certified_gap_bound(0, 1e-3, 2).bound > 0      # a zero defect is allowed
        # A seed outside the admissible region is not certified.
        @test certified_gap_bound(1 // 4, 1 // 4, 1).bound === nothing
    end

    @testset "rational chain with the bootstrap tail" begin
        record = rational_gap_chain(frac(1, 100), frac(1, 100), 2304, frac(784); steps = 6)
        claim = parse_rational(record["gap_strict_lower"])
        @test 0 < claim < certified_gap_bound(1 // 100, 1 // 100, 784).bound
        @test claim > log(20) / 784
        @test length(record["finite_intervals"]) == 6
        # Gaps above 1 are certified (the running minimum once started at 1).
        big_gap = rational_gap_chain(frac(1, 100), frac(1, 100), 2, frac(1); steps = 4)
        @test parse_rational(big_gap["gap_strict_lower"]) > 3
        @test bootstrap_coefficient(frac(0)) == frac(9, 2)
        w = beta_window(1e-4, 0.4)
        @test w["window_nonempty"] && w["beta_min_from_thermal_purity"] ≈ log(99) / 0.4
    end

    @testset "exact arithmetic" begin
        for x in (frac(1), frac(3, 2), frac(10), frac(12345, 7))
            lo, hi = log_enclosure(x)
            @test lo <= big(log(BigFloat(x))) <= hi && hi - lo < frac(1, 10^20)
        end
        lo, hi = power_interval(frac(2), frac(1, 2))
        @test lo^2 <= 2 <= hi^2 && hi - lo <= frac(1, 10^24)
        @test integer_root(big(10)^30 + 1, 3) == big(10)^10
        @test round_up(frac(1, 3), 100) == frac(34, 100)
        @test round_down(frac(-1, 3), 100) == frac(-34, 100)
        a, b = first_true(x -> x^2 >= 2, frac(1), frac(2); digits = 6)
        @test b - a == frac(1, 10^6) && a^2 < 2 <= b^2
        for x in (frac(3, 1024), frac(-5, 3), frac(7), big(1) // big(2)^200)
            @test parse_fraction(encode_fraction(x)) == x
        end
        @test encode_fraction(big(1) // big(2)^200) == "1//2^200"
        @test parse_rational("7/8") == parse_rational("7//8") == frac(7, 8)
        @test parse_rational(Dict("numerator" => "3", "denominator" => 4)) == frac(3, 4)
        @test_throws ErrorException parse_rational(0.5)
    end

    @testset "results and provenance" begin
        mktempdir() do dir
            path = joinpath(dir, "sub", "out.json")
            @test write_result(path, Dict("n" => 1)) == path
            @test endswith(read(path, String), "\n") && !isfile(path * ".tmp")
            # JSON3 would read a 40-digit integer as Float64; the exact parser keeps every digit.
            @test parse_exact_json("{\"n\": $(big(10)^40 + 1)}")["n"] == string(big(10)^40 + 1)
            mkpath(joinpath(dir, "src"))
            write(joinpath(dir, "src", "a.jl"), "x = 1\n")
            write(joinpath(dir, "Project.toml"), "")
            @test collect(keys(source_hashes(dir))) == [joinpath("src", "a.jl")]
            hashes = source_hashes(dir; files = ("Project.toml",))
            @test sort(collect(keys(hashes))) == ["Project.toml", joinpath("src", "a.jl")]
            @test hashes[joinpath("src", "a.jl")] == sha256_file(joinpath(dir, "src", "a.jl"))
            meta = run_metadata(Dict("k" => 1); source_sha256 = hashes, status = "test")
            @test meta["certified"] == false && meta["config"]["k"] == 1
        end
    end

    @testset "spin operators" begin
        for s in (1 // 2, 1, 3 // 2, 2)
            S = spin_matrices(s)
            @test S.sx * S.sy - S.sy * S.sx ≈ im * S.sz
            @test S.sx^2 + real(S.sy^2) + S.sz^2 ≈ s * (s + 1) * I
            @test S.sz[1, 1] == s                              # basis m = s, …, −s
        end
        plus, minus = integer_ladders(2)
        w = binomial_weights(2)
        @test w == [1, 4, 6, 4, 1]
        @test plus * Diagonal(w) == Diagonal(w) * transpose(minus)
        for s in (1, 2)
            G = cartesian_generators(s)
            @test all(g -> g ≈ -transpose(g), G)
            m = heisenberg(s)
            @test eigvals(Symmetric(m.rp_bond)) ≈ eigvals(Symmetric(m.chain.bonds(2)[1][3]))
            λ, anti = rp_test(m.rp_bond, m.d)
            @test λ > -1e-12 && anti < 1e-12
        end
        m = heisenberg(1)
        @test first(rp_test(m.chain.bonds(2)[1][3], 3)) ≈ -2       # standard basis is not RP
        Cs = rp_decomposition(m.rp_bond, m.d)
        @test norm(sum(kron(C, conj(C)) for C in Cs) + m.rp_bond) < 1e-12
    end

    @testset "exact diagonalisation" begin
        m = heisenberg(1)
        s = spectrum(m.chain, 6)
        @test s.exact
        @test s.E[1] ≈ -8.617423181814 atol = 1e-10
        rows = []
        flipped = spectrum(m.chain, 6; spin_flip = true, on_sector = r -> push!(rows, r))
        @test flipped.E ≈ s.E && flipped.N == 3^6
        @test [r.charge for r in rows] == 0:2:12 && all(r.multiplicity == 2 for r in rows[2:end])
        open = Chain(3, open_bonds(m.chain.bonds(2)[1][3]); charge = [2, 0, -2])
        @test spectrum(open, 2).E ≈ [-2; fill(-1, 3); fill(1, 5)]
        triangle = Chain(2, graph_bonds([(1, 2), (2, 3), (1, 3)], -kron(PAULI.z, PAULI.z)))
        @test spectrum(triangle, 3).E ≈ [-3, -3, 1, 1, 1, 1, 1, 1]
        truncated = PartialSpectrum(s.E[1:50], s.N, false)
        lo, hi = logZ_bounds(truncated, 0.5)
        exact, _ = logZ_bounds(s, 0.5)
        @test lo <= exact <= hi && exact ≈ log_partition(s.E, 0.5)
        stag = staggered_heisenberg(1.0)
        rotated = Chain(2, ring_bonds(stag.rp_bond); onsite = -spin_matrices(1 // 2).sz)
        @test spectrum(stag.chain, 6).E ≈ spectrum(rotated, 6).E atol = 1e-12
    end

    @testset "TFIM: exact diagonalisation against free fermions" begin
        t = tfim(2.0)
        for L in (4, 6, 8)
            ff = setprecision(BigFloat, 256) do
                Float64(tfim_logZ(L, 1.3, 2.0))
            end
            @test logZ_bounds(spectrum(t.chain, L), 1.3)[1] ≈ ff rtol = 1e-12
        end
        εₛ, εₜ = map(ε -> Float64(ε, RoundUp), tfim_seed_defects(8, 3, 2.0))
        @test certified_gap_bound(εₛ, εₜ, 3).bound ≈ 1.138848 atol = 1e-6
        spectra = Dict(L => spectrum(t.chain, L) for L in (4, 8))
        @test best_seed(scan_seeds(spectral_defects(spectra), [4], [1.0, 2.0])) !== nothing
    end

    @testset "purification extension" begin
        n, β = 4, 0.5
        sites = siteinds("S=1", n)
        pairs = [(j, mod1(j + 1, n)) for j in 1:n]
        bond(i, j) = op("Sz", sites, i) * op("Sz", sites, j) +
                     (op("S+", sites, i) * op("S-", sites, j) +
                      op("S-", sites, i) * op("S+", sites, j)) / 2
        halfgate(sites, i, j, τ) = exp(-τ * bond(i, j))
        curve = purified_log_partition(sites, pairs, halfgate, [0, β]; dbeta = 0.025,
            maxdim = 81, cutoff = 1e-14, extra = logz -> Dict("R" => expm1(logz)))
        @test curve[0.0]["logZ"] ≈ n * log(3)
        exact = log_partition(spectrum(heisenberg(1).chain, n).E, β)
        @test curve[β]["logZ"] ≈ exact atol = 1e-4
        @test haskey(curve[β], "R") && !haskey(curve[0.0], "R")
    end
end

include(joinpath(@__DIR__, "ring_certificates.jl"))
include(joinpath(@__DIR__, "certified_mpo", "runtests.jl"))
