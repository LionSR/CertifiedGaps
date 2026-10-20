# Verify every certificate in this directory; any mismatch throws or fails a test.
#
#   OPENBLAS_NUM_THREADS=1 julia --project=julia models/spin1_largeD/check.jl
#   OPENBLAS_NUM_THREADS=1 julia --project=julia models/spin1_largeD/check.jl --full
#
# Stage 1: each results/certificate_D*.json is checked by GapCertificates.check_point against its
# point (model, RP data, seed, params, stage-1 source hashes, the Lean seed hypotheses) and its
# exact chain is replayed. The fast check then recomputes log Z and the seed defects from the
# stored, hash-pinned verified spectra that stage 2 reads, and rebuilds the L = 6 spectra from the
# model. Stage 2: TransferCertificate.check_certificate replays results/certificate.json from the
# same spectra and recomputes the L = 8, 10 bridge. `--full` re-verifies the L = 6 and L = 12
# spectra of both points from scratch (check_point with spectra = true).
using Test
using CertifiedGaps.Results: read_exact_json
using CertifiedGaps.GapCertificates: check_point, seed_bounds, lean_points
using CertifiedGaps.VerifiedSpectra: verified_ring_spectrum
isdefined(@__MODULE__, :TransferCertificate) ||
    include(joinpath(@__DIR__, "src", "TransferCertificate.jl"))
using .TransferCertificate.LargeDModel
const TC = TransferCertificate
const FULL = "--full" in ARGS

@testset "Stage 1: core certificates at n0 = 6" begin
    lean = Dict(p.params[1] => p for p in lean_points(LEAN_FILE))
    @test sort(collect(keys(lean))) == sort([p.D for p in POINTS])
    for (D, beta0) in POINTS
        model = ring_model(D)
        @test lean[D].params == [D, beta0]
        path = joinpath(TC.ROOT, certificate_file(D))
        seconds = @elapsed recomputed = check_point(path, model; n0 = N0, β0 = beta0,
            params = Dict("D" => D), model_dir = TC.ROOT, sources = STAGE1_SOURCES,
            rp = rp_data(model), lean = lean[D], spectra = FULL)
        cert = read_exact_json(path)
        @test recomputed["gap_lower_bound"] == cert["gap_lower_bound"]
        # Stage 2 reads the stored spectra (hash-checked by reference_spectrum), so they must
        # reproduce the certificate's log Z enclosures and defects exactly.
        a, b = TC.reference_spectrum(D, N0), TC.reference_spectrum(D, 2N0)
        stored = seed_bounds(model, N0, beta0; spectra = Dict(N0 => a.spectrum,
            2N0 => b.spectrum))
        @test stored.log_partition == cert["log_partition"]
        @test (string(stored.eps_s), string(stored.eps_t)) ==
              (cert["seed"]["eps_s"], cert["seed"]["eps_t"])
        println("D=$D (n0, beta0) = ($N0, $beta0): eps_s <= $(Float64(stored.eps_s)), ",
            "eps_t <= $(Float64(stored.eps_t)), gap > $(cert["gap_lower_bound_decimal"])",
            FULL ? ", full replay in $(round(seconds; digits = 1)) s" : "")
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
