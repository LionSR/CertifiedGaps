# Produces the dimerized Heisenberg certificates in results/ with the shared core:
#
#   julia certify.jl rp          exact RP audit at the six deltas             -> rp_audit.json
#   julia certify.jl delta p q   one dimerization δ = p/q: verified spectra of the rings of 2, 4
#                                and 8 cells, a core gap certificate per seed of SEEDS
#                                (certificate_delta_*.json), the all-β exclusions at the other
#                                n₀ ∈ {2, 4}, the 16-spin gap, and its row of summary.json
#   julia certify.jl             all of the above, every δ
#
# The certificates are `GapCertificates.certify_point` outputs for `rotated_model(δ)`;
# check.jl verifies everything written here.
using CertifiedGaps, CertifiedGaps.Results, IntervalArithmetic, LinearAlgebra
using CertifiedGaps.RingModels, CertifiedGaps.VerifiedSpectra, CertifiedGaps.GapCertificates
using CertifiedGaps.Exclusions: exclude_all_beta
using CertifiedGaps.ExactArithmetic: round_up, round_down, parse_rational
isdefined(@__MODULE__, :DimerizedHeisenberg) ||
    include(joinpath(@__DIR__, "src", "DimerizedHeisenberg.jl"))
isdefined(@__MODULE__, :DimerizedRP) || include(joinpath(@__DIR__, "src", "rp.jl"))
isdefined(@__MODULE__, :ReflectionConventions) ||
    include(joinpath(@__DIR__, "src", "reflection_conventions.jl"))
using .DimerizedHeisenberg, .DimerizedRP
BLAS.set_num_threads(1)

const RESULTS_DIR = joinpath(@__DIR__, "results")
const SUMMARY = joinpath(RESULTS_DIR, "summary.json")
const GRID = big(10)^6                  # outward decimal grid of the reported gap intervals

# Run metadata with the model-local and core SHA-256 of a non-certificate result.
metadata(config; certified) = model_metadata(MODEL_DIR, config;
    packages = (CertifiedGaps, IntervalArithmetic), certified)

_slug(x) = denominator(x) == 1 ? "$(numerator(x))" : "$(numerator(x))_$(denominator(x))"

"Path of the core certificate of the seed (δ, β₀)."
certificate_path(δ, β0) =
    joinpath(RESULTS_DIR, "certificate_delta_$(_slug(δ))_beta_$(_slug(β0)).json")

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
            "real_identical_factor_identity" => a.real_identical_identity,
            "onsite_real_symmetric" => true,
            "certified_gap_bound" => nothing, "rational_seeds" => nothing,
        ))
    end
    data = Dict(
        "metadata" => metadata(Dict("model" => "dimerized_heisenberg", "sites_per_cell" => 2);
            certified = false),
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
            "bond. The reflected form with the cell swap holds exactly (the certificates).",
        "largest_exact_audit_matrix_dimension" => 256)
    path = write_result(joinpath(RESULTS_DIR, "rp_audit.json"), data)
    println("Wrote $path; exact audit seconds = $elapsed")
end

"""
    certify_delta(δ) -> Dict

Everything for one dimerization from one verified spectrum per ring: the core certificate of
every seed in `SEEDS`, the all-β exclusions (`exclude_all_beta`) at the even n₀ ∈ {2, 4}
without a seed, the enclosed 16-spin gap, and the dimer perturbation bound. Writes the
certificates and returns the summary row.
"""
function certify_delta(δ)
    start = time()
    model, rp = rotated_model(δ), rotated_rp(δ)
    spectra = Dict(L => verified_ring_spectrum(model, L) for L in (2, 4, 8))
    rows = Any[]
    for (n0, β0) in get(SEEDS, δ, Tuple{Int,Rational{BigInt}}[])
        seed = seed_bounds(model, n0, β0; spectra = Dict(n0 => spectra[n0], 2n0 => spectra[2n0]))
        path = certificate_path(δ, β0)
        cert = certify_point(path, model; n0, β0, params = point_params(δ),
            model_dir = MODEL_DIR, packages = (IntervalArithmetic,), rp, seed)
        row = certificate_row(cert)
        row["file"], row["sha256"] = "results/" * basename(path), sha256_file(path)
        push!(rows, row)
        println("delta=$δ n0=$n0 beta0=$β0 eps_s=$(Float64(seed.eps_s)) ",
            "eps_t=$(Float64(seed.eps_t)) gap>=$(cert["interval_bound_decimal"]) ",
            "chain>$(cert["gap_lower_bound_decimal"])")
    end
    exclusions = [exclude_all_beta(spectra[n], spectra[2n], n) for n in (2, 4)
                  if !any(r -> r["n0"] == n, rows)]
    lo, hi = spectral_gap_enclosure(spectra[8])
    lo > 0 || error("the 16-spin ground state is not separated")
    gap = [round_down(lo, GRID), round_up(hi, GRID)]
    # The headline bound is the interval bound (`certified_gap_bound`), as in the root README.
    best = isempty(rows) ? nothing :
        maximum(Rational{BigInt}(r["interval_bound_decimal"]) for r in rows)
    ratio = best === nothing ? nothing :
        string.([round_down(best / gap[2], GRID), round_up(best / gap[1], GRID)])
    dimer = dimer_bounds(δ, 4, DIMER_BETA(δ))
    endpoints(x) = string.([Rational{BigInt}(inf(x)), Rational{BigInt}(sup(x))])
    return Dict("delta" => string(δ), "certificates" => rows, "negative_results" => exclusions,
        "reference_gap_16_interval" => string.(gap), "bound_to_gap16_interval" => ratio,
        "dimer_perturbation" => Dict("n0" => "4", "beta0" => string(DIMER_BETA(δ)),
            "norm_small" => string(dimer.norm_small), "eps_s_enclosure" => endpoints(dimer.s),
            "eps_t_enclosure" => endpoints(dimer.t)),
        "spectra" => Dict(string(L) => spectrum_summary(s) for (L, s) in spectra),
        "seconds" => time() - start)
end

"""
    certify_summary(deltas)

Regenerate `deltas` with `certify_delta` and write results/summary.json, keeping the stored rows
of the other dimerizations.
"""
function certify_summary(deltas)
    old = isfile(SUMMARY) ? read_result(SUMMARY)["parameters"] : Any[]
    rows = Dict(parse_rational(r["delta"]) => r for r in old)
    for δ in deltas
        rows[δ] = certify_delta(δ)
    end
    Set(keys(rows)) == Set(DELTAS) || error("summary rows missing; run certify.jl")
    certified = [δ for δ in DELTAS if !isempty(rows[δ]["certificates"])]
    write_result(SUMMARY, Dict(
        "format" => "dimerized-heisenberg-summary-v4", "parameters" => [rows[δ] for δ in DELTAS],
        "model" => "H = sum over a ring of 2L spins of (1+delta) S.S inside each two-spin " *
            "cell plus (1-delta) S.S between neighboring cells",
        "frame" => "Certificates are for the rotated frame U'HU, U = (1 x J)^{x L}, " *
            "J = [0 1; -1 0] on the second spin of every cell (rotated_model, frame_identity)",
        "smallest_certified_grid_delta" => string(minimum(certified)),
        "validity" => "every even L >= n0 = 4 cells (N = 2L spins divisible by 4, N >= 8)",
        "negative_scope" => "For delta=1/5,3/20,1/10 and n0=2,4 cells (n0=2 at every delta), " *
            "every beta>0 fails the doubling recurrence of certified_gap_bound. This does not " *
            "imply gaplessness or exclude other methods, other seed theorems or larger rings.",
        "ratio_scope" => "Certified lower bound divided by the enclosed 16-spin gap; not an " *
            "infinite-chain ratio",
        "metadata" => metadata(Dict("model" => "dimerized_heisenberg"); certified = true)))
    println("Wrote $SUMMARY")
end

if abspath(PROGRAM_FILE) == @__FILE__
    step = isempty(ARGS) ? "all" : ARGS[1]
    step in ("all", "rp", "delta") || error("usage: certify.jl [rp | delta p q]")
    step in ("all", "rp") && certify_rp()
    if step == "delta"
        length(ARGS) == 3 || error("usage: certify.jl delta p q")
        certify_summary([parse(BigInt, ARGS[2]) // parse(BigInt, ARGS[3])])
    elseif step == "all"
        certify_summary(DELTAS)
    end
end
