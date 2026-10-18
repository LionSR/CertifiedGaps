# Fast tests: model frame, reflection and RP, verified spectra against independent ED, the fast
# checks of all 30 certificates, tamper detection, one full replay and the seed table.
using Test, LinearAlgebra
include(joinpath(@__DIR__, "..", "src", "LongitudinalCertificate.jl"))
using .LongitudinalCertificate, CertifiedGaps, CertifiedGaps.Results, CertifiedGaps.RingModels
using CertifiedGaps.VerifiedSpectra: verified_ring_spectrum, log_partition_enclosure
using CertifiedGaps.ExactDiagonalization: Chain, ring_bonds, spectrum
using CertifiedGaps.Models: tfim_logZ
using IntervalArithmetic: in_interval
const C = LongitudinalCertificate

# Independently written Pauli matrices.
const PX, PZ = Rational{BigInt}[0 1; 1 0], Rational{BigInt}[1 0; 0 -1]
levels(spec) = sort(reduce(vcat, [repeat(s.centers, s.multiplicity) for s in spec]))

@testset "Frame, reflection and exact RP" begin
    for (g, h) in C.points()
        m = C.longitudinal(g, h)
        # The certificate frame is the physical one: the frame unitary is the identity.
        @test is_frame_change(RingModel(2, -kron(PZ, PZ), -g * PX - h * PZ), m, I(2))
        @test m.reflection == I(2)
        @test verify_rp(m, C.rp(m)) && verify_rp(m, rp_certificate(m))
        @test isreal_model(m) && all(iszero, m.charge)
    end
    # The site reflection must map the bond to its mirror image: a non-symmetric bond fails.
    @test_throws ArgumentError RingModel(2, -kron(PZ, PX) - kron(PX, PZ) + kron(PZ, PX);
        reflection = I(2))
end

@testset "Verified spectra against computational-basis ED" begin
    for (g, h) in ((1//4, 1//20), (1//1, 1//4), (2//1, 1//1)), L in (4, 6)
        spec = verified_ring_spectrum(C.longitudinal(g, h), L)
        chain = Chain(2, ring_bonds(-Float64.(kron(PZ, PZ)));
            onsite = Float64.(-g * PX - h * PZ))
        @test maximum(abs, levels(spec) - spectrum(chain, L; dense_max = 2^L).E) < 1e-12
    end
    # h = 0: the enclosure contains the free-fermion log Z.
    z = log_partition_enclosure(verified_ring_spectrum(C.longitudinal(2//1, 0//1), 8), 3//2)
    @test in_interval(Float64(tfim_logZ(8, 3//2, 2)), z)
end

lean = C.lean_seeds()
@testset "Fast checks of all certificates" begin
    @test C.check_grid(C.result_paths())
    for (g, h) in C.points()
        @test C.check(g, h; lean)["gap_lower_bound"] ==
              read_result(C.certificate_path(g, h))["gap_lower_bound"]
    end
end

@testset "Tampering is detected" begin
    g, h = 1//1, 1//20
    dir = mktempdir()
    # A copy of the stored certificate with one field changed, at a path of the same name.
    function tampered(f)
        c = read_result(C.certificate_path(g, h))
        f(c)
        path = joinpath(dir, basename(C.certificate_path(g, h)))
        write_result(path, c)
        path
    end
    fails(f; full = false, lean = lean) =
        (@test_throws Exception C.check(g, h; path = tampered(f), full, lean))
    @test C.check(g, h; path = tampered(identity), lean)["gap_lower_bound"] ==
          read_result(C.certificate_path(g, h))["gap_lower_bound"]
    fails(c -> c["seed"]["eps_s"] = "1//100")
    fails(c -> c["gap_lower_bound"] = "1//10")
    fails(c -> c["interval_bound_decimal"] = 0.08)
    fails(c -> c["seed"]["beta0"] = "7//1")
    fails(c -> c["params"]["h"] = "1//10")
    fails(c -> c["model"]["onsite"]["re"][1][1] = "-1//10")
    fails(c -> delete!(c["model"], "reflection"))
    fails(c -> c["metadata"]["source_sha256"]["models/tfim_longitudinal/check.jl"] = "0"^64)
    # Stored defects above the Lean seed defects fail even when the certificate is consistent.
    fails(identity; lean = Dict((g, h) => (s = 18//100, t = 1//100)))
    # A stored εₛ below the fresh enclosure (≈ 0.18399987) fails the full replay.
    fails(c -> c["seed"]["eps_s"] = "18399//100000"; full = true)
end

@testset "Full replay of the critical point g = 1, h = 1/20" begin
    r = C.check(1//1, 1//20; full = true, lean)
    @test r["defect_s_upper"] <= 0.184 && r["defect_t_upper"] <= 0.0075
end

@testset "Seed table" begin
    @test C.check_seed_table()
    @test length(read_result(joinpath(C.ROOT, "results", "seeds.json"))["points"]) == 30
    @test length(lean) == 30
end
