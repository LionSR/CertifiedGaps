using CertifiedGaps, LinearAlgebra, Test
using CertifiedGaps.Results
isdefined(@__MODULE__, :check_historical_hashes) || include("provenance.jl")
BLAS.set_num_threads(1)
include("rp.jl")
using .DimerizedRP
function exact_rational(s)
    parts = split(s, "//")
    return parse(BigInt, parts[1]) // parse(BigInt, parts[2])
end
function check_audit()
    data = read_exact_json(joinpath(@__DIR__, "results", "rp_audit.json"))
    @test data["exact_algebraic_audit"] === true
    @test data["gap_certificate"] === false
    @test data["metadata"]["certified"] === false
    @test length(data["parameters"]) == 6
    for row in data["parameters"]
        a = audit_rp(exact_rational(row["delta"]))
        @test a.defect_norm2 == exact_rational(row["swap_defect_frobenius_squared"])
        @test a.real_rhs_error_norm2 == exact_rational(row["proposed_real_rhs_error_frobenius_squared"])
        @test a.defect_norm2 > 0
        @test row["conjugate_identity_exact"] === true
        @test a.reflected_real_identity === true
        @test row["reflected_real_identity_exact"] === true
        @test row["theta_invariant_onsite_exact"] === true
        @test row["repository_general_RP_eligible"] === true
        @test row["real_identical_factor_identity"] === false
        # There are no valid seed inputs; never manufacture or recompute a gap from them.
        @test row["rational_seeds"] === nothing
        @test row["certified_gap_bound"] === nothing
    end
    check_source_hashes(data["metadata"])
end
if abspath(PROGRAM_FILE) == @__FILE__
    @testset "Exact RP obstruction audit (not a gap certificate)" begin
        check_audit()
    end
end
