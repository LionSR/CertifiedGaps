# Verify every certificate in this directory; any mismatch throws or fails a test.
#
#   OPENBLAS_NUM_THREADS=1 julia --project=models/spin1_largeD models/spin1_largeD/check.jl
#   OPENBLAS_NUM_THREADS=1 julia --project=models/spin1_largeD models/spin1_largeD/check.jl --full
#
# Stage 1 (test/ed_seeds.jl): the ED seeds results/certificate_D*.json from their stored,
# hash-checked spectral enclosures. Stage 2 (TransferCertificate.check_certificate): every case
# of results/certificate.json, plus a fresh recomputation of the L = 8, 10 bridge spectra.
# The fast check also rebuilds the two L = 6 enclosures from the model (under a second each).
# `--full` adds the dense L = 12 replays (about 35 min per D, one BLAS thread).
using Test
include("test/ed_seeds.jl")
isdefined(@__MODULE__, :TransferCertificate) || include("src/TransferCertificate.jl")

@testset "Transfer certificate and finite bridge" begin
    @test TransferCertificate.check_certificate() !== nothing
end

"""
    replay_enclosure(D, L)

Rebuild the complete verified spectrum at (D, L) with `compute_spectrum` from `model.jl` and
compare it with the stored, hash-checked one. Both enclose the same sorted spectrum, each
sorted centre within its list's largest radius of the true level, so the centres must agree to
within the sum of the two largest radii; this holds whether or not LAPACK reproduces the stored
bits. Returns the number of bit-identical centres.
"""
function replay_enclosure(D, L)
    stored = TransferCertificate.reference_spectrum(D, L)
    fresh = TransferCertificate.compute_spectrum(L, D; rigorous = true)
    @test fresh["complete"] && fresh["rigorous"]
    E, R = Float64.(fresh["energies"]), Float64.(fresh["radii"])
    @test length(E) == length(stored.E) == 3^L
    tolerance = maximum(R) + maximum(stored.R)
    @test all(abs(a - b) <= tolerance for (a, b) in zip(E, stored.E))
    same = count(E .== stored.E)
    println("replay D=$D L=$L: $same of $(3^L) centres bit-identical, max radius ",
        maximum(R), " (stored ", maximum(stored.R), ")")
    return same
end

@testset "Fresh dense replay of the L = 6 enclosures" begin
    for D in (4//1, 1603//500)
        replay_enclosure(D, 6)
    end
end

if "--full" in ARGS
    @testset "Fresh dense replay of the L = 12 enclosures" begin
        for D in (4//1, 1603//500)
            replay_enclosure(D, 12)
        end
    end
end
