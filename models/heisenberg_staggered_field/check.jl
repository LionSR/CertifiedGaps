# Checker for the staggered-field Heisenberg certificates.
#
#     julia --project=models/heisenberg_staggered_field \
#         models/heisenberg_staggered_field/check.jl            # fast, ~45 s
#     ... check.jl --full                                        # + ~6 min
#
# The fast check (test/runtests.jl) recomputes log Z, the rational seeds and the bounds of all
# five certificates from the stored eigenvalues and radii, and the exact parts of the
# shared-core cross-check. `--full` also re-verifies the cross-check's spectra from the model
# with CertifiedGaps.GapCertificates. scripts/spectral_replay.jl rebuilds every block of the
# five main certificates (about 3-4 min) and rewrites results/spectral_replay_summary.json.
include("test/runtests.jl")
if "--full" in ARGS
    using CertifiedGaps.GapCertificates: check_certificate
    @testset "shared-core cross-check at h = 1 (verified spectra)" begin
        cert = CertifiedGaps.Results.read_exact_json(
            joinpath(ROOT, "results", "core_certificate_1_1.json"))
        recomputed = check_certificate(cert; spectra = true, progress = true)
        @test recomputed["gap_lower_bound"] == cert["gap_lower_bound"]
    end
end
