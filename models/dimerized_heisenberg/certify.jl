# Produces the dimerised Heisenberg certificate files in results/:
#
#   julia certify.jl rp           exact RP audit at the six deltas        -> rp_audit.json
#   julia certify.jl seeds        seed certificate or all-beta exclusion at every delta
#                                 (about 1 min each)                       -> seed_delta_*.json.gz
#   julia certify.jl seeds 3 10   one delta only
#   julia certify.jl summary      plain index of the six certificates     -> certificate.json
#   julia certify.jl              all of the above, in this order
#
# check.jl verifies everything written here. A rerun records hashes of the current sources;
# the eigenvector solves need not be bitwise reproducible, so a rerun may store slightly
# different (equally valid) enclosures. Seed upper bounds are rounded up to a 1e-6 grid plus one
# spare unit, so that the stored seeds survive a fresh replay on another BLAS (check.jl --full).
using CertifiedGaps, CertifiedGaps.Results, IntervalArithmetic, LinearAlgebra
isdefined(@__MODULE__, :check_source_hashes) || include(joinpath(@__DIR__, "src", "provenance.jl"))
isdefined(@__MODULE__, :ReflectionConventions) ||
    include(joinpath(@__DIR__, "src", "reflection_conventions.jl"))
isdefined(@__MODULE__, :Certification) ||
    include(joinpath(@__DIR__, "src", "Certification.jl"))
using .Certification, .Certification.DimerizedRP
BLAS.set_num_threads(1)

const RESULTS_DIR = joinpath(@__DIR__, "results")
const MODEL = "models/dimerized_heisenberg"

"Exact RP audit of all six deltas; writes results/rp_audit.json (no gap claim)."
function certify_rp()
    rows = Any[]
    elapsed = @elapsed for delta in DELTAS
        a = audit_rp(delta)
        push!(rows, Dict(
            "delta" => string(delta), "weak_bond" => string(a.weak),
            "strong_bond" => string(a.strong),
            "swap_defect_frobenius_squared" => string(a.defect_norm2),
            "proposed_real_rhs_error_frobenius_squared" => string(a.real_rhs_error_norm2),
            "conjugate_identity_exact" => a.conjugate_identity,
            "reflected_real_identity_exact" => a.reflected_real_identity,
            "theta_invariant_onsite_exact" => a.theta_invariant_onsite,
            "repository_general_RP_eligible" => true,
            "real_identical_factor_identity" => a.real_identical_identity,
            "onsite_real_symmetric" => true,
            "certified_gap_bound" => nothing, "rational_seeds" => nothing,
            "seed_scan" => "see results/seed_delta_*.json.gz",
        ))
    end
    # Model-local sources are enforced by check.jl; the core is recorded as provenance only.
    hashes = source_hashes(repository_root(); dirs=(), files=(
        "$MODEL/src/rp.jl", "$MODEL/src/reflection_conventions.jl", "$MODEL/certify.jl",
        "$MODEL/check.jl", "$MODEL/Project.toml", "$MODEL/Manifest.toml"))
    data = Dict(
        "metadata" => merge(run_metadata(Dict("model" => "dimerized_heisenberg",
                "sites_per_cell" => 2); source_sha256=hashes, packages=(CertifiedGaps,),
            status="exact_obstruction_to_requested_real_RP_form", certified=false),
            Dict("core_sha256" => core_provenance())),
        "exact_algebraic_audit" => true, "gap_certificate" => false,
        "arithmetic" => "Complex{Rational{BigInt}}; normalized basis evaluated with exact " *
            "denominators",
        "audit_seconds" => elapsed, "parameters" => rows,
        "alternative_conventions" => [Dict("block_spins" => r.block_spins,
            "local_dimension" => r.local_dimension,
            "unit_boundary_swap_defect" => string(r.unit_boundary_swap_defect),
            "reflected_form" => r.reflected_form,
            "uniform_identical_form" => r.uniform_identical_form)
            for r in ReflectionConventions.audit_conventions()],
        "scope" => "No common cell basis admits sum c*C tensor C for the stated inter-cell " *
            "bond. The broader conjugate form is verified exactly and is not excluded.",
        "largest_exact_audit_matrix_dimension" => 256,
        "validity_range" => nothing, "reference_gap" => nothing,
        "bound_to_gap_ratio" => nothing)
    path = write_result(joinpath(RESULTS_DIR, "rp_audit.json"), data)
    println("Wrote $path; exact audit seconds = $elapsed")
end

"""
    certify_summary()

Write results/certificate.json: for every delta the SHA-256 of its seed certificate, the seeds
(or the all-beta exclusions), the enclosed 16-spin gap and the bound-to-gap16 ratio interval,
and the smallest certified grid delta. check.jl compares every entry with the certificates.
"""
function certify_summary()
    rows = Any[]
    for delta in DELTAS
        path = seed_path(delta)
        record = read_result(path)
        lo, hi = parse_rational.(record["reference_gap_16_interval"])
        gap = interval(inf(CertifiedGaps.to_interval(lo)), sup(CertifiedGaps.to_interval(hi)))
        best = isempty(record["seeds"]) ? nothing :
            maximum(parse_rational(s["gap_lower"]) for s in record["seeds"])
        ratio = best === nothing ? nothing :
            endpoints(CertifiedGaps.to_interval(best) / gap)
        push!(rows, Dict("delta" => string(delta), "source_file" => "results/" * basename(path),
            "sha256" => sha256_file(path), "seeds" => record["seeds"],
            "survey_seed" => record["survey_seed"], "seconds_total" => record["seconds_total"],
            "reference_gap_16_interval" => record["reference_gap_16_interval"],
            "bound_to_gap16_interval" => ratio,
            "negative_results" => [Dict("n0_cells" => r["n0_cells"], "status" => r["status"],
                "beta_boxes" => string(length(r["cover"]))) for r in record["negative_results"]]))
        println("delta=$delta bound=$(best === nothing ? "none" : Float64(best)) ratio=$ratio")
    end
    certified = [parse_rational(r["delta"]) for r in rows if !isempty(r["seeds"])]
    write_result(joinpath(RESULTS_DIR, "certificate.json"), Dict(
        "format" => "dimerized-heisenberg-summary-v2", "parameters" => rows,
        "model" => "H = sum over a ring of 2n spins of (1+delta) S.S inside each two-spin " *
            "cell plus (1-delta) S.S between neighbouring cells",
        "smallest_certified_grid_delta" => string(minimum(certified)),
        "maximum_certification_sites" => "16",
        "validity_sites" => "N divisible by 4, N >= 8",
        "negative_scope" => "For delta=1/5,3/20,1/10 and n0=2,4 cells, every beta>0 fails the " *
            "doubling recurrence of certified_gap_bound. This does not imply gaplessness or " *
            "exclude other methods, other seed theorems or larger rings.",
        "ratio_scope" => "Certified lower bound divided by the enclosed 16-spin gap; not an " *
            "infinite-chain ratio",
        "metadata" => merge(run_metadata(Dict("model" => "dimerized_heisenberg");
            source_sha256=source_hashes(repository_root(); dirs=(),
                files=("$MODEL/certify.jl", "$MODEL/src/Certification.jl")),
            packages=(CertifiedGaps, IntervalArithmetic),
            status="summary_of_seed_certificates", certified=true),
            Dict("core_sha256" => core_provenance()))))
    println("Wrote results/certificate.json")
end

if abspath(PROGRAM_FILE) == @__FILE__
    step = isempty(ARGS) ? "all" : ARGS[1]
    step in ("all", "rp", "seeds", "summary") || error("unknown step $step")
    step in ("all", "rp") && certify_rp()
    if step == "seeds" && length(ARGS) == 3
        generate(parse(BigInt, ARGS[2]) // parse(BigInt, ARGS[3]))
    elseif step in ("all", "seeds")
        foreach(generate, DELTAS)
    end
    step in ("all", "summary") && certify_summary()
end
