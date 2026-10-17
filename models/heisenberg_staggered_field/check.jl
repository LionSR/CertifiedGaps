# Checker for the staggered-field Heisenberg certificates (core format, one per field).
#
#     julia --project=julia models/heisenberg_staggered_field/check.jl   # exact parts
#     julia --project=julia models/heisenberg_staggered_field/check.jl --full   # + spectra
#
# The fast check runs test/runtests.jl: exact frame, RP data, blocking, the exact chain of every
# certificate (`check_certificate(cert; spectra = false)`), the Lean seed constants, source
# hashes and the summary. `--full` reruns `check_certificate(cert; spectra = true)`, which
# rebuilds and re-verifies the 8- and 16-site spectra of every point and requires the fresh
# defects to be at most the stored εₛ, εₜ.
include(joinpath(@__DIR__, "test", "runtests.jl"))
if "--full" in ARGS
    @testset "verified seed spectra (full replay)" begin
        for h in FIELDS
            cert = read_exact_json(certificate_path(h))
            seconds = @elapsed recomputed = check_certificate(cert; spectra = true)
            @test recomputed["gap_lower_bound"] == cert["gap_lower_bound"]
            @test recomputed["defect_s_upper"] <= parse_rational(cert["seed"]["eps_s"])
            @test recomputed["defect_t_upper"] <= parse_rational(cert["seed"]["eps_t"])
            println("h = $h: spectra re-verified in $(round(seconds; digits = 1)) s")
        end
    end
end
