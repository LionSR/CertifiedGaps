# Produces the dimerised Heisenberg certificate files in results/:
#
#   julia certify.jl rp          exact RP audit at six deltas         -> rp_audit.json
#   julia certify.jl verified    verified seeds at delta = 1/2, 3/10, 1/4
#                                (about a minute each)               -> verified_delta_*.json.gz
#   julia certify.jl verified 3 10    one delta only
#   julia certify.jl summary     plain copy of the three seeds        -> certificate.json
#   julia certify.jl             all of the above, in this order
#
# check.jl verifies everything written here. A rerun records hashes of the current sources;
# the eigenvector solves need not be bitwise reproducible (see scripts/independent_review.jl).
using CertifiedGaps, CertifiedGaps.Results, IntervalArithmetic, JSON3, LinearAlgebra
include(joinpath(@__DIR__, "src", "provenance.jl"))
include(joinpath(@__DIR__, "src", "rp.jl"))
include(joinpath(@__DIR__, "src", "reflection_conventions.jl"))
include(joinpath(@__DIR__, "src", "verified_spectra.jl"))
using .DimerizedRP, .VerifiedSpectra
BLAS.set_num_threads(1)

const RESULTS_DIR = joinpath(@__DIR__, "results")
const MODEL = "models/dimerized_heisenberg"
"The deltas of the RP audit (the three certified ones and three smaller ones)."
const AUDIT_DELTAS = [big(1)//2, big(3)//10, big(1)//4, big(1)//5, big(3)//20, big(1)//10]
"The certified deltas and their seed temperatures β₀ (from the numerical scan)."
const SEEDS = Dict(big(1)//2 => big(6)//1, big(3)//10 => big(13)//2, big(1)//4 => big(29)//4)

"Exact RP audit of all six deltas; writes results/rp_audit.json (no gap claim)."
function certify_rp()
    rows = Any[]
    elapsed = @elapsed for delta in AUDIT_DELTAS
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
            "seed_scan" => "see results/verified_delta_*.json.gz and numerical_summary.json",
        ))
    end
    hashes = source_hashes(repository_root(); dirs=("julia/src",), files=(
        "$MODEL/src/rp.jl", "$MODEL/src/reflection_conventions.jl", "$MODEL/certify.jl",
        "$MODEL/check.jl", "julia/Project.toml", "julia/Manifest.toml",
        "outputs/thermal-mpo-gap-scope.tex"))
    data = Dict(
        "metadata" => run_metadata(Dict("model" => "dimerized_heisenberg",
                "sites_per_cell" => 2); source_sha256=hashes, packages=(CertifiedGaps,),
            status="exact_obstruction_to_requested_real_RP_form", certified=false),
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
    certify_verified(delta, beta)

Enclose the complete spectra of the 8- and 16-spin rings, form the outward log-partition
enclosures, the exact rational seeds εₛ ≥ 1 − Pₛ(4, β), εₜ ≥ 1 − P_τ(8, β) and the uniform
bound; writes results/verified_delta_<p>_<q>.json.gz.
"""
function certify_verified(delta::Rational{BigInt}, beta::Rational{BigInt})
    0 <= delta < 1 && beta > 0 || error("need 0 <= delta < 1 and beta > 0")
    audit_rp(delta).reflected_real_identity || error("reflected RP identity fails")
    t = time()
    small, Esmall, rsmall = verified_ring(delta, 8)
    large, Elarge, rlarge = verified_ring(delta, 16)
    lsmall = enclosed_logZ(Esmall, rsmall, beta)
    llarge = enclosed_logZ(Elarge, rlarge, beta)
    ldouble = enclosed_logZ(Elarge, rlarge, 2beta)
    epss = rational_upper(1 - exp(llarge - 2lsmall))
    epst = rational_upper(1 - exp(ldouble - 2llarge))
    0 <= epss < 1 && 0 <= epst < 1 || error("seed defects outside [0, 1)")
    bound = certified_gap_bound(epss, epst, beta)
    bound.bound === nothing && error("seed not admissible at delta = $delta, beta = $beta")
    sources = source_hashes(repository_root(); dirs=("julia/src",), files=(
        "$MODEL/src/rp.jl", "$MODEL/src/sectors.jl", "$MODEL/src/verified_spectra.jl",
        "$MODEL/certify.jl", "$MODEL/check.jl", "$MODEL/test/verified_independent.jl",
        "$MODEL/METHOD.md", "julia/Project.toml", "julia/Manifest.toml",
        "outputs/thermal-mpo-gap-note.tex", "outputs/thermal-mpo-gap-scope.tex"))
    record = Dict("metadata" => run_metadata(Dict("delta" => string(delta),
                "beta0" => string(beta), "n0_cells" => 4);
            source_sha256=sources, packages=(CertifiedGaps, IntervalArithmetic),
            status="verified_spectral_seeds_general_reflected_RP", certified=true),
        "certification_scope" => "General reflected RP theorem Theta=F_cell K; NOT the " *
            "original stronger real identical-factor representation",
        "original_stronger_RP_condition_satisfied" => false,
        "general_reflected_RP_exact" => true, "delta" => string(delta), "n0_cells" => 4,
        "beta0" => string(beta), "epsilon_s_upper" => string(epss),
        "epsilon_t_upper" => string(epst),
        "log_partition_enclosures" => Dict(
            "n_beta" => string.([rational_lower(lsmall), rational_upper(lsmall)]),
            "2n_beta" => string.([rational_lower(llarge), rational_upper(llarge)]),
            "2n_2beta" => string.([rational_lower(ldouble), rational_upper(ldouble)])),
        "recurrence_bound" => bound.bound,
        "recurrence_bound_exact_binary" => string(Rational{BigInt}(bound.bound)),
        "recurrence_diagnostics" => bound,
        "validity_range" => "even L >= 4 cells; 8,12,16,... physical spins",
        "seconds_total" => time() - t, "rings" => [small, large])
    path = joinpath(RESULTS_DIR,
        "verified_delta_$(numerator(delta))_$(denominator(delta)).json")
    write_result(path, record)
    # Stored compressed (about 0.5 MB instead of 2.7 MB); check.jl reads the .gz.
    run(`gzip -9 -n -f $path`)
    println("SAVED $path.gz epsilon_s=$(Float64(epss)) epsilon_t=$(Float64(epst)) ",
        "bound=$(bound.bound) seconds=$(time() - t)")
end

"Seed data of every stored verified certificate, as plain exact rationals."
function certificate_summary()
    paths = sort(filter(p -> startswith(basename(p), "verified_delta_") &&
        endswith(p, ".json.gz"), readdir(RESULTS_DIR; join=true)))
    rows = Any[]
    for path in paths
        data = JSON3.read(read_json_text(path))
        push!(rows, Dict(
            "delta" => String(data.delta), "n0_cells" => Int(data.n0_cells),
            "n0_spins" => 2 * Int(data.n0_cells), "beta0" => String(data.beta0),
            "epsilon_s_upper" => String(data.epsilon_s_upper),
            "epsilon_t_upper" => String(data.epsilon_t_upper),
            "log_partition_enclosures" => Dict(String(k) => String.(v)
                for (k, v) in pairs(data.log_partition_enclosures)),
            "recurrence_bound_exact_binary" => String(data.recurrence_bound_exact_binary),
            "validity_range" => String(data.validity_range),
            "source_file" => "results/" * basename(path),
            "source_file_sha256" => sha256_file(path)))
    end
    return rows
end

"Write results/certificate.json, the seeds of the three verified certificates without spectra."
function certify_summary()
    hashes = source_hashes(repository_root(); dirs=(),
        files=("$MODEL/certify.jl", "$MODEL/src/provenance.jl"))
    path = write_result(joinpath(RESULTS_DIR, "certificate.json"), Dict(
        "metadata" => run_metadata(Dict("model" => "dimerized_heisenberg",
                "sites_per_cell" => 2); source_sha256=hashes, packages=(CertifiedGaps,),
            status="summary_of_verified_certificates", certified=true),
        "model" => "H = sum over a ring of 2n spins of (1+delta) S.S inside each two-spin " *
            "cell plus (1-delta) S.S between neighbouring cells",
        "criterion_lengths" => "n counts two-spin cells; P_s(n0,beta0) = Z_{2n0}/Z_{n0}^2, " *
            "P_tau(2n0,beta0) = Z_{2n0}(2beta0)/Z_{2n0}(beta0)^2",
        "rp_form" => "Scope-note bond form with an internal reflection: rotate the second " *
            "spin of every cell by J = [0 1; -1 0]; then -h = (1-delta) sum_{C in (Sx, iSy, " *
            "Sz) on the second spin of the left cell} C (x) Theta(C), Theta(X) = F conj(X) F " *
            "with F the swap of the two spins of a cell, and o = (1+delta) J-rotated S1.S2 " *
            "is real with Theta(o) = o. In the cell basis M/sqrt(2) (F = M M^T/2, columns " *
            "|00>+|11>, |00>-|11>, |01>+|10>, i(|01>-|10>)) this is -h = sum C' (x) " *
            "conj(C') with complex C' and real diagonal o. No common cell basis gives a " *
            "swap-symmetric sum C (x) C. Exact checks: src/rp.jl.",
        "certificates" => certificate_summary()))
    println("Wrote $path")
end

if abspath(PROGRAM_FILE) == @__FILE__
    step = isempty(ARGS) ? "all" : ARGS[1]
    step in ("all", "rp") && certify_rp()
    if step == "verified" && length(ARGS) == 3
        delta = parse(BigInt, ARGS[2]) // parse(BigInt, ARGS[3])
        certify_verified(delta, SEEDS[delta])
    elseif step in ("all", "verified")
        for delta in sort(collect(keys(SEEDS)); rev=true)
            certify_verified(delta, SEEDS[delta])
        end
    end
    step in ("all", "summary") && certify_summary()
    step in ("all", "rp", "verified", "summary") || error("unknown step $step")
end
