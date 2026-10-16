# Checker for the dimerised Heisenberg certificates.
#
# Fast (about 1 min): from the stored files it
#   1. re-derives every exact RP identity of results/rp_audit.json at all six deltas;
#   2. for each results/seed_delta_*.json.gz: checks the source hashes, rebuilds the model and
#      verifies the stored RP factors exactly, re-derives every block's residual, orthogonality
#      and radius witness from its stored component maxima, checks sector completeness and the
#      momentum-conjugate blocks, recomputes every seed's three log-partition enclosures, defects
#      and `certified_gap_bound`, the 16-spin gap interval, the dimer perturbation bound, and
#      every all-beta exclusion (both analytic tails and each interval box of the cover);
#   3. checks that results/certificate.json indexes those certificates exactly.
# `check.jl --full` (about 4 min) then re-solves every block of every ring from the model,
# rebuilds the eigenvectors, residuals and Gram defects, and repeats all of step 2 with the
# fresh spectra only: the stored eigenvalues and residual maxima are not trusted.
using CertifiedGaps, CertifiedGaps.Results, IntervalArithmetic, LinearAlgebra, Test
isdefined(@__MODULE__, :check_source_hashes) || include(joinpath(@__DIR__, "src", "provenance.jl"))
isdefined(@__MODULE__, :Certification) ||
    include(joinpath(@__DIR__, "src", "Certification.jl"))
using .Certification, .Certification.VerifiedSpectra, .Certification.DimerizedRP
using .Certification: exact # IntervalArithmetic also exports `exact`.
BLAS.set_num_threads(1)

const RESULTS_DIR = joinpath(@__DIR__, "results")

"Parse an exact `\"num//den\"` string."

"Re-derive the exact RP audit (reflected form holds, identical-factor form fails)."
function check_audit()
    data = read_exact_json(joinpath(RESULTS_DIR, "rp_audit.json"))
    @test data["exact_algebraic_audit"] === true
    # The audit is about the RP hypothesis only; the seeds live in the seed certificates.
    @test data["gap_certificate"] === false
    @test data["metadata"]["certified"] === false
    @test [parse_rational(row["delta"]) for row in data["parameters"]] == DELTAS
    for row in data["parameters"]
        a = audit_rp(parse_rational(row["delta"]))
        @test a.defect_norm2 == parse_rational(row["swap_defect_frobenius_squared"])
        @test a.real_rhs_error_norm2 ==
              parse_rational(row["proposed_real_rhs_error_frobenius_squared"])
        @test a.defect_norm2 > 0
        @test row["conjugate_identity_exact"] === true
        @test a.reflected_real_identity === true
        @test row["reflected_real_identity_exact"] === true
        @test row["theta_invariant_onsite_exact"] === true
        @test row["repository_general_RP_eligible"] === true
        @test row["real_identical_factor_identity"] === false
        @test row["rational_seeds"] === nothing
        @test row["certified_gap_bound"] === nothing
    end
    check_source_hashes(data["metadata"])
end

"""
    check_seed_certificates(; fresh=false)

Check every seed certificate with `Certification.check_record`: from the stored spectra after
re-deriving each block witness (`fresh=false`), or from spectra re-solved from the model
(`fresh=true`, check.jl --full).
"""
function check_seed_certificates(; fresh=false)
    for delta in DELTAS
        record = read_result(seed_path(delta))
        @test record["delta"] == string(delta)
        spectra = fresh ? nothing : stored_spectra(record)
        seconds = @elapsed @test check_record(record; spectra) isa Dict
        bounds = [Float64(parse_rational(s["gap_lower"])) for s in record["seeds"]]
        println("checked delta=$delta fresh=$fresh bounds=$bounds ",
            "exclusions=$(length(record["negative_results"])) seconds=$(round(seconds; digits=1))")
    end
end

"results/certificate.json must index the six seed certificates exactly."
function check_certificate_summary()
    summary = read_result(joinpath(RESULTS_DIR, "certificate.json"))
    check_source_hashes(summary["metadata"])
    @test summary["format"] == "dimerized-heisenberg-summary-v2"
    @test [parse_rational(r["delta"]) for r in summary["parameters"]] == DELTAS
    positive = Rational{BigInt}[]
    for row in summary["parameters"]
        path = seed_path(parse_rational(row["delta"]))
        @test row["source_file"] == "results/" * basename(path)
        @test sha256_file(path) == row["sha256"]
        record = read_result(path)
        for key in ("delta", "seeds", "survey_seed", "reference_gap_16_interval")
            @test record[key] == row[key]
        end
        @test [(r["n0_cells"], r["status"], string(length(r["cover"])))
            for r in record["negative_results"]] ==
            [(r["n0_cells"], r["status"], r["beta_boxes"]) for r in row["negative_results"]]
        if isempty(record["seeds"])
            @test row["bound_to_gap16_interval"] === nothing
        else
            push!(positive, parse_rational(row["delta"]))
            bound = maximum(parse_rational(s["gap_lower"]) for s in record["seeds"])
            gl, gh = parse_rational.(record["reference_gap_16_interval"])
            rl, rh = parse_rational.(row["bound_to_gap16_interval"])
            @test rl <= bound / gh <= bound / gl <= rh
        end
    end
    # Every smaller grid delta has no seed; check_record proved its n0 = 2, 4 exclusions.
    @test minimum(positive) == parse_rational(summary["smallest_certified_grid_delta"])
    @test all(d -> d in positive, filter(>=(minimum(positive)), DELTAS))
end

"""
    refined_witness(rec, L)

Rigorous residual, orthogonality and eigenvalue-radius bounds for one freshly solved block
record `rec` of an `L`-spin ring, sharper than `VerifiedSpectra.scalar_witness` (Frobenius
factor 3d/2 instead of 2d). Used by scripts/independent_review.jl.
"""
function refined_witness(rec, L)
    d = Int(rec["dimension"])
    h = parse_rational(rec["matrix_error_upper"])
    q = parse_rational(rec["q_component_max"])
    lam = parse_rational(rec["lambda_max_abs"])
    rawr = parse_rational(rec["raw_residual_component_max"])
    rawe = parse_rational(rec["raw_orthogonality_component_max"])
    coarse = VerifiedSpectra.scalar_witness(d, L, parse_rational(rec["h_component_max"]),
        q, lam, rawr, rawe, h)
    @assert coarse.ebound < 1 # Thus ||Q|| < 2, independently of the refinement.
    a = big(3) * L // 4
    # Cauchy on two concatenated real dot-product vectors gives these sums.
    rerr = min(coarse.rerr, coarse.gamma * (2 * (a + h) + q * lam) + coarse.underflow)
    eerr = min(coarse.eerr, coarse.gamma * 5 + coarse.underflow)
    # ||complex dxd matrix||F <= sqrt(2)*d*componentmax < (3/2)*d*componentmax.
    factor = big(3) * d // 2
    r = factor * (rawr + rerr); e = factor * (rawe + eerr)
    @assert e < 1
    return (; r, e, radius = r + 2 * h + (a + lam) * e)
end

if abspath(PROGRAM_FILE) == @__FILE__
    @testset "Dimerised Heisenberg certificate" begin
        @testset "Exact RP audit" begin
            check_audit()
        end
        @testset "Seed certificates and exclusions (stored spectra)" begin
            check_seed_certificates()
        end
        @testset "Summary index" begin
            check_certificate_summary()
        end
        if "--full" in ARGS
            @testset "Fresh spectral replay of every certificate" begin
                check_seed_certificates(; fresh=true)
            end
        end
    end
end
