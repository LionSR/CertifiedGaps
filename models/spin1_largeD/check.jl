# Verify every certificate in this directory; any mismatch throws or fails a test.
#
#   OPENBLAS_NUM_THREADS=1 julia --project=models/spin1_largeD models/spin1_largeD/check.jl
#
# Stage 1 (test/ed_seeds.jl): the ED seeds results/certificate_D*.json from their stored,
# hash-checked spectral enclosures. Stage 2 (TransferCertificate.check_certificate): every case
# of certificate.json, plus a fresh recomputation of the L = 8, 10 bridge spectra.
# Neither stage repeats the 35-minute L = 12 diagonalisations (see README.md).
using Test
include("test/ed_seeds.jl")
isdefined(@__MODULE__, :TransferCertificate) || include("src/TransferCertificate.jl")

@testset "Transfer certificate and finite bridge" begin
    @test TransferCertificate.check_certificate() !== nothing
end
