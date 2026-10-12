# This produces an exact obstruction audit, NOT a gap certificate.
using ThermalMPOGap, LinearAlgebra
using ThermalMPOGap.Results
BLAS.set_num_threads(1)
include("rp.jl")
using .DimerizedRP
include("reflection_conventions.jl")
const ROOT = normpath(joinpath(@__DIR__, "../.."))
const DELTAS = [big(1)//2, big(3)//10, big(1)//4, big(1)//5, big(3)//20, big(1)//10]
function main()
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
            "seed_scan" => "see_separate_numerical_delta_results",
        ))
    end
    hashes = source_hashes(ROOT; dirs=("julia/src", "certificates/dimerized_heisenberg/test"),
        files=("certificates/dimerized_heisenberg/rp.jl",
            "certificates/dimerized_heisenberg/certify.jl",
            "certificates/dimerized_heisenberg/check.jl",
            "certificates/dimerized_heisenberg/reflection_conventions.jl",
            "certificates/dimerized_heisenberg/numerical_ed.jl",
            "julia/Project.toml", "julia/Manifest.toml",
            "outputs/thermal-mpo-gap-scope.tex"))
    data = Dict(
        "metadata" => run_metadata(Dict("model"=>"dimerized_heisenberg", "sites_per_cell"=>2);
            source_sha256=hashes, packages=(ThermalMPOGap,),
            status="exact_obstruction_to_requested_real_RP_form", certified=false),
        "exact_algebraic_audit" => true, "gap_certificate" => false,
        "arithmetic" => "Complex{Rational{BigInt}}; normalized basis evaluated with exact denominators",
        "audit_seconds" => elapsed, "parameters" => rows,
        "alternative_conventions" => [Dict("block_spins"=>r.block_spins,
            "local_dimension"=>r.local_dimension,
            "unit_boundary_swap_defect"=>string(r.unit_boundary_swap_defect),
            "reflected_form"=>r.reflected_form,
            "uniform_identical_form"=>r.uniform_identical_form)
            for r in ReflectionConventions.audit_conventions()],
        "scope" => "No common cell basis admits sum c*C tensor C for the stated inter-cell bond. The broader conjugate form is verified exactly and is not excluded.",
        "largest_exact_audit_matrix_dimension" => 256,
        "numerical_spectra" => "see numerical_delta_*.json; not certified",
        "validity_range" => nothing, "reference_gap" => nothing,
        "bound_to_gap_ratio" => nothing,
        "finite_size_seed_extrapolation" => "see numerical_summary.json for unvalidated broader-criterion seed extrapolation",
    )
    path = write_result(joinpath(@__DIR__, "results", "rp_audit.json"), data)
    println("Wrote $path; exact audit seconds = $elapsed; no gap certificate")
end
main()
