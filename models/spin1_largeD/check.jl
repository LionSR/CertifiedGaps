# Verify every certificate in this directory; any mismatch throws or fails a test.
#
#   OPENBLAS_NUM_THREADS=1 julia --project=julia models/spin1_largeD/check.jl
#   OPENBLAS_NUM_THREADS=1 julia --project=julia models/spin1_largeD/check.jl --full
#
# Stage 1: each results/certificate_D*.json is a core certificate. The fast check verifies its
# model-local source hashes, that its model, RP data and seed (n0, beta0) are the ones defined in
# src/model.jl, its exact RP identity and rational chain (GapCertificates.check_certificate with
# spectra = false), and recomputes log Z and the seed defects from the stored, hash-pinned
# verified spectra; it also rebuilds the L = 6 spectra from the model. Stage 2:
# TransferCertificate.check_certificate replays results/certificate.json from the same spectra
# and recomputes the L = 8, 10 bridge. `--full` re-verifies the L = 6 and L = 12 spectra of both
# points from scratch (GapCertificates.check_certificate with spectra = true).
using Test, IntervalArithmetic
using CertifiedGaps.Results: read_exact_json, verify_model_metadata
using CertifiedGaps.RingModels: decode_model, decode_rp, encode_rp
using CertifiedGaps.VerifiedSpectra: verified_ring_spectrum, log_partition_enclosure
import CertifiedGaps.GapCertificates
isdefined(@__MODULE__, :TransferCertificate) ||
    include(joinpath(@__DIR__, "src", "TransferCertificate.jl"))
using .TransferCertificate.LargeDModel
const TC = TransferCertificate

"Exact equality of two `RingModel`s, field by field."
same_model(a, b) = all(getfield(a, f) == getfield(b, f) for f in fieldnames(typeof(a)))

@testset "Stage 1: core certificates at n0 = 6" begin
    for (D, beta0) in POINTS
        cert = read_exact_json(joinpath(TC.ROOT, certificate_file(D)))
        @test verify_model_metadata(cert["metadata"], TC.ROOT; STAGE1_SOURCES...) == 2
        model = ring_model(D)
        @test same_model(decode_model(cert["model"]), model)
        @test encode_rp(decode_rp(cert["reflection_positivity"])) == encode_rp(rp_data(model))
        seed = cert["seed"]
        Q = Rational{BigInt}
        @test parse(Int, seed["n0"]) == N0 && parse(Q, seed["beta0"]) == beta0
        es, et = parse(Q, seed["eps_s"]), parse(Q, seed["eps_t"])
        @test GapCertificates.check_certificate(cert; spectra = false)["gap_lower_bound"] ==
              cert["gap_lower_bound"]
        # The hypotheses of lean/CertifiedGaps/Models/SpinOneLargeD.lean.
        @test es <= LEAN_SEED_BOUNDS[D][1] && et <= LEAN_SEED_BOUNDS[D][2]
        # log Z and the defects from the stored spectra (hash-checked by reference_spectrum).
        a, b = TC.reference_spectrum(D, N0), TC.reference_spectrum(D, 2N0)
        zn, z2n = TC.logz(a, beta0), TC.logz(b, beta0)
        z2n2 = TC.logz(b, 2beta0)
        stored = cert["log_partition"]
        for (key, z) in (("logZ_n0_beta0", zn), ("logZ_2n0_beta0", z2n),
            ("logZ_2n0_2beta0", z2n2))
            @test parse.(Q, stored[key]) == [Q(inf(z)), Q(sup(z))]
        end
        @test sup(1 - exp(z2n - 2zn)) <= es          # 1 - P_s(n0, beta0)
        @test sup(1 - exp(z2n2 - 2z2n)) <= et        # 1 - P_tau(2 n0, beta0)
        println("D=$D (n0, beta0) = ($N0, $beta0): eps_s <= $(Float64(es)), ",
            "eps_t <= $(Float64(et)), gap > $(cert["gap_lower_bound_decimal"])")
    end
end

@testset "Fresh L = 6 spectra against the stored ones" begin
    for (D, _) in POINTS
        stored = TC.reference_spectrum(D, N0).spectrum
        fresh = verified_ring_spectrum(ring_model(D), N0)
        @test length(fresh) == length(stored)
        for (f, s) in zip(fresh, stored)
            @test f.multiplicity == s.multiplicity && length(f.centers) == length(s.centers)
            # Both enclose the same sorted block spectrum, whatever bits LAPACK returns.
            @test all(abs(x - y) <= f.radius + s.radius for (x, y) in zip(f.centers, s.centers))
        end
    end
end

@testset "Stage 2: transfer certificate and finite bridge" begin
    @test TC.check_certificate() !== nothing
end

if "--full" in ARGS
    @testset "Full replay: verified L = 6 and L = 12 spectra from the model" begin
        for (D, _) in POINTS
            cert = read_exact_json(joinpath(TC.ROOT, certificate_file(D)))
            seconds = @elapsed recomputed =
                GapCertificates.check_certificate(cert; spectra = true, progress = true)
            @test recomputed["gap_lower_bound"] == cert["gap_lower_bound"]
            println("D=$D: full replay in $(round(seconds; digits = 1)) s")
        end
    end
end
