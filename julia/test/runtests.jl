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
        # The other reference values formerly printed by julia/scripts/constants.jl.
        B(x) = interval(BigFloat, x)
        setprecision(BigFloat, 256) do
            # Theorem 2.22 (i) at ε₀ = 1/100: β₀Δ_L ≥ log(0.97/(3ε₀)) − Λ₄₀ ≈ 3.3471.
            @test 3.34707 < inf(log(B(97) / 3) - Λ) < 3.34708
            # Theorem 2.22 (ii), r₀ ≤ 0.03 and a₀ ≤ 1/100.
            for (R₀, low, high) in ((B(3) / 100, 0.833876, 0.833877),
                                    (B(1) / 99, 0.787726, 0.787727))
                Λii, okii = majorant_loss(R₀, B(1); J = 40)
                @test okii && low < sup(Λii) < high
            end
        end
        # Odd lengths (Proposition 2.13): β₀Δ_L ≥ log((1 − 2u₀)/(2u₀)), u₀ = 1 − 0.99³.
        u₀ = 1 - (B(99) / 100)^3
        @test 2.762187 < inf(log((1 - 2u₀) / (2u₀))) < 2.762188
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

    @testset "compressed and compact results, model provenance" begin
        mktempdir() do dir
            data = Dict("n" => typemax(Int64), "x" => "3//4", "v" => [1, 2])
            pretty = write_result(joinpath(dir, "a.json"), data)
            compact = write_result(joinpath(dir, "b.json"), data; compact = true)
            packed = write_result(joinpath(dir, "c.json.gz"), data)
            @test read(compact, String) == CertifiedGaps.Results.JSON3.write(data) * "\n"
            @test read_result_bytes(packed) == read(pretty)          # gzip of the pretty text
            @test read(`gzip -dc -- $packed`) == read(pretty) && !isfile(packed * ".tmp")
            @test read_result(packed) == read_result(pretty) == read_exact_json(pretty)
            @test read_result(packed)["n"] == string(typemax(Int64))         # exact token
            @test read_result(packed; exact = false)["v"] == [1, 2]
            @test read_result(packed, Dict{String,Any})["x"] == "3//4"
            # Equal JSON gives equal compressed bytes (gzip -n: no name, no timestamp).
            @test read(write_result(joinpath(dir, "d.json.gz"), data)) == read(packed)
            @test_throws ArgumentError write_result(joinpath(dir, "e.json"), data; gzip = true)
            @test_throws ArgumentError write_result(joinpath(dir, "e.gz"), data; gzip = false)

            repo = joinpath(dir, "repo")
            model = joinpath(repo, "models", "toy")
            for (path, text) in (("models/toy/src/a.jl", "a = 1\n"),
                ("models/toy/check.jl", "c = 1\n"), ("models/toy/Project.toml", ""),
                ("models/toy/Manifest.toml", ""), ("julia/src/Core.jl", "k = 1\n"),
                ("julia/Project.toml", ""), ("julia/Manifest.toml", ""))
                mkpath(dirname(joinpath(repo, path)))
                write(joinpath(repo, path), text)
            end
            p = model_provenance(model; files = ("check.jl", "Project.toml", "Manifest.toml"))
            @test sort(collect(keys(p.source_sha256))) == ["models/toy/Manifest.toml",
                "models/toy/Project.toml", "models/toy/check.jl", "models/toy/src/a.jl"]
            @test sort(collect(keys(p.core_sha256))) ==
                  ["julia/Manifest.toml", "julia/Project.toml", "julia/src/Core.jl"]
            recorded = merge(p.source_sha256, p.core_sha256,
                Dict("../outside.jl" => "0"^64))
            @test verify_source_hashes(recorded, repo; scope = "models/toy/") == 4
            # Editing the core does not matter in the model scope, editing the model does.
            write(joinpath(repo, "julia/src/Core.jl"), "k = 2\n")
            @test verify_source_hashes(recorded, repo; scope = "models/toy/") == 4
            @test_throws ErrorException verify_source_hashes(recorded, repo)
            write(joinpath(repo, "models/toy/src/a.jl"), "a = 2\n")
            @test_throws ErrorException verify_source_hashes(recorded, repo;
                scope = "models/toy/")
            @test_throws ErrorException verify_source_hashes(recorded, repo; scope = "none/")
            # Model-relative keys work with the model directory as root.
            local_keys = source_hashes(model; dirs = ("src",))
            @test verify_source_hashes(local_keys, model) == 1
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
        # Exact spin-½ and Pauli matrices.
        σ, S = exact_pauli(), exact_spin_half()
        @test σ.x * σ.x == σ.id && σ.iy == im * σ.y && σ.x * σ.z == -σ.iy
        @test eltype(σ.x) == Rational{BigInt} && eltype(exact_pauli(Int).z) == Int
        @test S.x * S.y - S.y * S.x == im * S.z
        @test S.x^2 + S.y^2 + S.z^2 == 3 // 4 * S.id
        @test kron(S.x, S.x) - kron(S.iy, S.iy) + kron(S.z, S.z) ≈
              real(sum(kron(A, A) for A in (spin_matrices(1 // 2).sx, spin_matrices(1 // 2).sy,
            spin_matrices(1 // 2).sz)))
        # Integer spin-1 Cartesian generators: S = −iG is a spin-1 representation.
        G = integer_cartesian_generators()
        @test eltype(G[1]) == Int && all(g -> transpose(g) == -g, G)
        @test Float64.(G[1]) == cartesian_generators(1)[1] && G == map(g -> Int.(g),
            cartesian_generators(1))
        Sg = [-im * g for g in G]
        @test Sg[1] * Sg[2] - Sg[2] * Sg[1] == im * Sg[3]
        @test sum(s * s for s in Sg) == 2I
        @test -sum(kron(s, s) for s in Sg) == sum(kron(g, g) for g in G)
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
        # Weighted spectra: multiplicities instead of repeated levels.
        weighted = WeightedSpectrum([0.0, 1.0, 2.0], [1, 3, 5], 9, true)
        expanded = PartialSpectrum([0.0; ones(3); fill(2.0, 5)], 9, true)
        for β in (0.125, 1.0, 8.0)
            wl, wh = logZ_bounds(weighted, β)
            @test wl == wh
            @test wl ≈ logZ_bounds(expanded, β)[1] atol = 1e-14
        end
        partial = WeightedSpectrum([0.0, 1.0], [1, 3], 9, false)
        pl, ph = logZ_bounds(partial, 1.0)
        @test pl < logZ_bounds(weighted, 1.0)[1] <= ph
        stag = staggered_heisenberg(1.0)
        rotated = Chain(2, ring_bonds(stag.rp_bond); onsite = -spin_matrices(1 // 2).sz)
        @test spectrum(stag.chain, 6).E ≈ spectrum(rotated, 6).E atol = 1e-12
    end

    @testset "rigorous seed defects from log Z enclosures" begin
        t = tfim(2.0)
        spectra = Dict(L => spectrum(t.chain, L) for L in (4, 8))
        εₛ, εₜ = seed_defects(spectra, 4, 1.5)
        # Thin intervals around the float values, widened by a margin for rounding.
        w(x) = interval(x - 1e-12, x + 1e-12)
        zn = w(logZ_bounds(spectra[4], 1.5)[1])
        z2n = w(logZ_bounds(spectra[8], 1.5)[1])
        z2n2 = w(logZ_bounds(spectra[8], 3.0)[1])
        r = seed_defects(zn, z2n, z2n2)
        @test r.eps_s isa Rational{BigInt} && r.eps_t isa Rational{BigInt}
        @test inf(r.defect_s) <= εₛ <= r.eps_s && inf(r.defect_t) <= εₜ <= r.eps_t
        @test r.eps_s <= Rational{BigInt}(sup(1 - exp(z2n - 2zn)))
        @test r.eps_t <= Rational{BigInt}(sup(1 - exp(z2n2 - 2z2n)))
        @test abs(Float64(r.eps_s) - εₛ) < 1e-9
        es, et = r                                   # positional destructuring: εₛ, εₜ
        @test (es, et) == (r.eps_s, r.eps_t)
        # BigFloat inputs: the bound equals the directed-rounding −expm1 of the lower end.
        setprecision(BigFloat, 256) do
            a, b = interval(BigFloat, 1 // 10), interval(BigFloat, 3 // 20)
            s = seed_defects(a, b, b)
            x = setrounding(BigFloat, RoundDown) do
                BigFloat(3 // 20, RoundDown) - 2 * BigFloat(1 // 10, RoundUp)
            end
            @test s.eps_s == Rational{BigInt}(-setrounding(() -> expm1(x), BigFloat, RoundDown))
        end
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
