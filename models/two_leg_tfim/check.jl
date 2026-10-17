# Checker for the two-leg TFIM certificates (results/certificate_g*.json, core format).
#
#     julia --project=julia models/two_leg_tfim/check.jl          # fast
#     julia --project=julia models/two_leg_tfim/check.jl --full   # + spectra
#
# The fast check (test/runtests.jl) re-verifies the exact RP identity, the frame change, the
# stored model against `Ladder.ladder(g)`, the exact rational chain from the stored seeds
# (`GapCertificates.check_certificate(cert; spectra = false)`) and the model-local source hashes.
# `--full` also recomputes every verified spectrum and the seed defects from the model
# (`check_certificate(cert; spectra = true)`), which must not exceed the stored εₛ, εₜ.
include(joinpath(@__DIR__, "test", "runtests.jl"))
if "--full" in ARGS
    @testset "Seed defects from freshly verified spectra" begin
        for path in certificate_paths()
            fresh = check_certificate(read_exact_json(path); spectra = true, progress = true)
            println(basename(path), ": defect upper ends ", fresh["defect_s_upper"], ", ",
                fresh["defect_t_upper"])
            @test fresh["gap_lower_bound"] == read_exact_json(path)["gap_lower_bound"]
        end
    end
end
