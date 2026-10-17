# Fast tests: model frame and RP, verified spectra against independent ED, the fast checks of
# all 30 certificates, tamper detection, one full replay, the seed table and the Lean seeds.
using Test, LinearAlgebra
include(joinpath(@__DIR__, "..", "src", "LongitudinalCertificate.jl"))
using .LongitudinalCertificate, CertifiedGaps, CertifiedGaps.Results, CertifiedGaps.RingModels
using CertifiedGaps.VerifiedSpectra: verified_ring_spectrum, log_partition_enclosure
using CertifiedGaps.ExactDiagonalization: Chain, ring_bonds, spectrum
using CertifiedGaps.ExactArithmetic: parse_rational
using CertifiedGaps.Models: tfim_logZ
using IntervalArithmetic: in_interval
const C = LongitudinalCertificate

# Independently written Pauli matrices.
const PX, PZ = Rational{BigInt}[0 1; 1 0], Rational{BigInt}[1 0; 0 -1]
levels(spec) = sort(reduce(vcat, [repeat(s.centers, s.multiplicity) for s in spec]))

@testset "Frame and exact RP" begin
    U = Matrix{Rational{BigInt}}(I, 2, 2)
    for (g, h) in C.points()
        m = C.longitudinal(g, h)
        # The certificate frame is the physical one: the frame unitary U is the identity, so
        # the stored terms must equal the conjugated physical bond and on-site term exactly.
        @test m.bond == kron(U, U)' * (-kron(PZ, PZ)) * kron(U, U)
        @test m.onsite == U' * (-g * PX - h * PZ) * U
        @test verify_rp(m, C.rp(m)) && verify_rp(m, rp_certificate(m))
        @test isreal_model(m) && all(iszero, m.charge)
    end
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

paths = C.result_paths()
@testset "Fast checks of all certificates" begin
    @test C.check_grid(paths)
    for path in paths
        @test C.check_point(path)["gap_lower_bound"] == read_result(path)["gap_lower_bound"]
    end
end

@testset "Tampering is detected" begin
    path = C.certificate_path(1//1, 1//20)
    name = basename(path)
    tampered(f) = (c = read_result(path); f(c); c)
    @test_throws Exception C.check_record(tampered(c -> c["seed"]["eps_s"] = "1//100"), name)
    @test_throws Exception C.check_record(tampered(c -> c["gap_lower_bound"] = "1//10"), name)
    @test_throws Exception C.check_record(tampered(c -> c["seed"]["beta0"] = "7//1"), name)
    @test_throws Exception C.check_record(tampered(c -> c["h"] = "1//10"), name)
    @test_throws Exception C.check_record(
        tampered(c -> c["model"]["onsite"]["re"][1][1] = "-1//10"), name)
    @test_throws Exception C.check_record(tampered(c ->
        c["metadata"]["source_sha256"]["models/tfim_longitudinal/check.jl"] = "0"^64), name)
    # A stored εₛ below the fresh enclosure (≈ 0.18399987) fails the full replay.
    low = tampered(c -> c["seed"]["eps_s"] = "18399//100000")
    @test_throws Exception C.check_record(low, name; full = true)
end

@testset "Full replay of the critical point g = 1, h = 1/20" begin
    r = C.check_point(C.certificate_path(1//1, 1//20); full = true)
    @test r["defect_s_upper"] <= 0.184 && r["defect_t_upper"] <= 0.0075
end

@testset "Seed table and the Lean seeds" begin
    table = read_result(joinpath(C.ROOT, "results", "seeds.json"))
    @test length(table["points"]) == 30
    @test table["points"] == parse_exact_json(Results.JSON3.write(C.seed_rows()))
    # lean/CertifiedGaps/Models/IsingLongitudinal.lean states each point's seed defects as the
    # first entries of its `s` and `t` lists (rounded up); the certified ones must not exceed
    # them, at the same (g, h, β₀).
    lean = read(joinpath(C.ROOT, "..", "..", "lean", "CertifiedGaps", "Models",
        "IsingLongitudinal.lean"), String)
    frac(s) = (p = split(s, " / "); parse(BigInt, p[1]) //
        (length(p) == 1 ? 1 : parse(BigInt, p[2])))
    sci(s) = (m = match(r"^(\d+)e-(\d+)$", s); parse(BigInt, m[1]) // big(10)^parse(Int, m[2]))
    rows = Dict((parse_rational(r["g"]), parse_rational(r["h"])) => r for r in table["points"])
    pattern = r"⟨([^,]+), ([^,]+), ([^,]+), \{\s*s := \[([^,]+),[^\]]*\]\s*t := \[([^,]+),"s
    found = 0
    for m in eachmatch(pattern, lean)
        r = rows[(frac(m[1]), frac(m[2]))]
        @test parse_rational(r["beta0"]) == frac(m[3]) && r["n0"] == "6"
        @test parse_rational(r["epsilon_s"]) <= sci(m[4])
        @test parse_rational(r["epsilon_t"]) <= sci(m[5])
        found += 1
    end
    @test found == 30
end
