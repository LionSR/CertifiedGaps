# Fast tests (about 1 min): the fast checker (check.jl), the exact frame and RP identities, the
# verified spectra of small rings against the closed-form 4-spin spectrum and the Float64 sector
# solver, the Lean seed constants, the dimer perturbation bound, the recurrence rejection and the
# checker's rejection of corrupted certificates.
include(joinpath(@__DIR__, "..", "check.jl"))
include(joinpath(@__DIR__, "..", "src", "reflection_conventions.jl"))
using CertifiedGaps.Recurrence: to_interval
const Q = Rational{BigInt}

@testset "Stored certificates (fast check)" begin
    @test check_audit()
    @test check_summary()
end

@testset "Exact rotated frame and reflected RP" begin
    for δ in DELTAS
        model = rotated_model(δ)
        # U = (1 ⊗ J)^{⊗L} maps the physical Hamiltonian onto the certified one, exactly.
        @test frame_identity(δ)
        @test FRAME' * FRAME == I
        @test verify_rp(model, rotated_rp(δ))
        # The LDL† factorization with the cell swap finds an RP certificate as well; without
        # the swap the bond is not reflection positive in this basis.
        @test verify_rp(model, rp_certificate(model; reflection = rotated_rp(δ).reflection))
        @test_throws ArgumentError rp_certificate(model)
        # A wrong frame (no rotation) does not give the certified bond.
        onsite, bond = physical_terms(δ)
        @test exact_matrix(bond) != model.bond
    end
end

@testset "Exact alternative cell conventions and the conjugated-factor form" begin
    conventions = ReflectionConventions.audit_conventions()
    @test length(conventions) == 2
    @test all(row.reflected_form && !row.uniform_identical_form for row in conventions)
    for δ in DELTAS
        a = audit_rp(δ)
        m = RingModel(4, a.transformed_bond, a.transformed_onsite)
        @test verify_rp(m, rp_certificate(m; factors = a.channels, weights = fill(1 - δ, 3)))
        @test a.basis_numerator' * a.basis_numerator == 2Matrix{Complex{Q}}(I, 4, 4)
    end
end

@testset "Verified four-spin spectrum and log Z against the closed form" begin
    for δ in DELTAS
        spectrum = verified_ring_spectrum(rotated_model(δ), 2)
        @test sum(length(s.centers) * s.multiplicity for s in spectrum) == 16
        # Two singlets, three triplets, one quintet; 256-bit intervals, independent of Float64.
        setprecision(BigFloat, 256) do
            x = interval(BigFloat, numerator(δ)) / interval(BigFloat, denominator(δ))
            root = sqrt(1 + 3x^2)
            exact = vcat([-1 - root, -1 + root], fill(interval(BigFloat, -1), 3), fill(-x, 3),
                fill(x, 3), fill(interval(BigFloat, 1), 5))
            sort!(exact; by = inf)
            balls = sort([(c, BigFloat(sup(to_interval(s.radius)))) for s in spectrum
                          for c in s.centers for _ in 1:s.multiplicity])
            for ((c, r), e) in zip(balls, exact)
                @test max(abs(inf(e) - BigFloat(c)), abs(sup(e) - BigFloat(c))) < r
            end
            z = log_partition_enclosure(spectrum, 6 // big(1))
            direct = log(sum(exp.(-6 .* exact)))
            @test BigFloat(inf(z)) <= inf(direct) <= sup(direct) <= BigFloat(sup(z))
        end
    end
end

include("numerical_blocks.jl")      # also defines the Float64 `dense_spectrum` of sectors.jl
@testset "Verified 8-spin spectra against the physical-frame Float64 solver" begin
    for δ in (1 // big(2), 1 // big(10))
        spectrum = verified_ring_spectrum(rotated_model(δ), 4)
        levels = sort(reduce(vcat, [repeat(s.centers, s.multiplicity) for s in spectrum]))
        reference, _ = dense_spectrum(Float64(δ), 8)
        @test maximum(abs, levels - reference.E) < 1e-12
    end
end

include(joinpath(@__DIR__, "..", "src", "Diagnostics.jl"))
@testset "Sparse momentum/parity blocks against complete spectra" begin
    for L in (4, 6, 8, 10)
        reference, _ = dense_spectrum(0.3, L)
        sparse, records = Diagnostics.low_spectrum(0.3, L; dense_max = 2000)
        @test sparse.exact
        @test maximum(abs, sparse.E - reference.E) < 1e-11
        @test sum(parse(Int, b["dimension"]) * parse(Int, b["multiplicity"])
                  for b in records) == 2^L
    end
end

@testset "Lean seed constants dominate the certified defects" begin
    # lean/CertifiedGaps/Models/DimerizedHeisenberg.lean takes 1 − Pₛ ≤ s[0], 1 − P_τ ≤ t[0]
    # at (n₀, β₀) = (4, β₀) as hypotheses; the certificates must supply them.
    lean = read(joinpath(REPO, "lean", "CertifiedGaps", "Models", "DimerizedHeisenberg.lean"),
        String)
    decimal(m, e) = parse(BigInt, m) // big(10)^parse(Int, e)
    fraction(s) = (p = parse.(BigInt, split(s, "/")); length(p) == 1 ? p[1] // 1 : p[1] // p[2])
    points = [(fraction(m[1]), fraction(replace(m[2], " " => "")), decimal(m[3], m[4]),
        decimal(m[5], m[6])) for m in eachmatch(
        r"⟨(\d+ / \d+), ([\d /]+), \{\s*s := \[(\d+)e-(\d+),[^⟩]*?t := \[(\d+)e-(\d+),", lean)]
    @test length(points) == 3
    summary = read_result(joinpath(RESULTS_DIR, "summary.json"))
    for (δ, β0, s0, t0) in points
        row = only(r for r in summary["parameters"] if parse_rational(r["delta"]) == δ)
        c = only(c for c in row["certificates"] if parse_rational(c["beta0"]) == β0)
        @test asint(c["n0"]) == 4
        @test parse_rational(c["eps_s"]) <= s0
        @test parse_rational(c["eps_t"]) <= t0
    end
end

@testset "Scalable dimer perturbation and recurrence rejection" begin
    d = dimer_bounds(999999 // big(1000000), 128, 8 // big(1))
    # A genuine beyond-ED example: largest ring 512 spins, no diagonalization.
    b = CertifiedGaps.certified_gap_bound(Q(sup(d.s)), Q(sup(d.t)), 8 // big(1)).bound
    @test b !== nothing && b > 0
    for δ in DELTAS
        bad = dimer_bounds(δ, 4, 6 // big(1))
        @test inf(bad.s) > 1 // 2
        @test rejects(Q(inf(bad.s)), Q(inf(bad.t))) !== nothing
    end
    @test rejects(1 // big(2), 0 // big(1)) == 0
    @test rejects(1 // big(100), 1 // big(100)) === nothing
end

@testset "Checker rejects corrupted records" begin
    summary = read_result(joinpath(RESULTS_DIR, "summary.json"))
    row = only(r for r in summary["parameters"] if r["delta"] == "1//2")
    @test check_delta(row)
    corrupt(f) = (bad = deepcopy(row); f(bad); bad)
    for f in (r -> r["certificates"][1]["eps_s"] = "0//1",
        r -> r["certificates"][1]["gap_lower_bound"] = "1//1",
        r -> r["certificates"][1]["n0"] = "2",
        r -> r["certificates"][1]["sha256"] = "0"^64,
        r -> popfirst!(r["negative_results"][1]["cover"]),
        r -> merge!(r["negative_results"][1]["cover"][1], Dict("s_lower" => "0//1",
            "t_lower" => "0//1")),
        r -> r["negative_results"][1]["tail_beta"] = "2//1",
        r -> r["negative_results"][1]["thermal_defect_lower_at_beta_1"] = "1//4",
        r -> empty!(r["negative_results"]),
        r -> r["bound_to_gap16_interval"][2] = "0//1",
        r -> r["dimer_perturbation"]["norm_small"] = "0//1")
        @test_throws ErrorException check_delta(corrupt(f))
    end
    # Tampering inside a certificate: its hash changes, and the core check fails on its own.
    cert = read_result(joinpath(@__DIR__, "..", row["certificates"][1]["file"]))
    @test check_certificate(cert; spectra = false) isa Dict
    for f in (c -> c["seed"]["eps_t"] = "0//1", c -> c["gap_lower_bound"] = "1//1",
        c -> c["reflection_positivity"]["factors"][1]["weight"] = "-1//1")
        bad = deepcopy(cert)
        f(bad)
        @test_throws Exception check_certificate(bad; spectra = false)
    end
    bad = deepcopy(summary)
    empty!(bad["metadata"]["source_sha256"])
    @test_throws ErrorException check_summary(; summary = bad)
end
