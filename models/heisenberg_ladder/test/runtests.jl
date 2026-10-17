# Fast tests: the exact frame and RP identities, core verified spectra against an independent
# dense Hamiltonian in the physical frame, the rejection witness, continuation, the all-β
# exclusion and its checker, and the stored suite and core certificates (no spectra).
using Test, LinearAlgebra, IntervalArithmetic
include(joinpath(@__DIR__, "..", "src", "HeisenbergLadder.jl"))
const HL = HeisenbergLadder
using CertifiedGaps.ExactDiagonalization: Chain, sector_hamiltonian
using CertifiedGaps.RingModels: RPCertificate, verify_rp, rp_certificate, encode_rp, decode_rp,
    exact_matrix
using CertifiedGaps.VerifiedSpectra: verified_ring_spectrum
using CertifiedGaps.Results: parse_exact_json, verify_model_metadata
using CertifiedGaps: to_interval

levels(spec) = sort(reduce(vcat, [repeat(s.centers, s.multiplicity) for s in spec]))

@testset "Exact frame unitary and reflection positivity" begin
    for J in HL.COUPLINGS
        @test HL.check_frame(J)
        model = HL.ladder_model(J)
        rp = HL.ladder_rp(model)
        @test verify_rp(model, rp) && rp.reflection == exact_matrix(HL.cell_reflection())
        bad = deepcopy(encode_rp(rp))
        bad["factors"][1]["weight"] = "2//1"
        @test !verify_rp(model, decode_rp(bad))
    end
    model = HL.ladder_model(4 // 1)
    # Without the cell reflection the bond is not RP in this frame (R[−h] is not PSD).
    @test_throws ArgumentError rp_certificate(model)
    # The exact LDL† factorisation with the reflection finds an equivalent certificate.
    @test verify_rp(model, rp_certificate(model; reflection = HL.cell_reflection()))
    G = exact_matrix(HL.leg_exchange())
    @test G * G == I && G == transpose(G)
end

@testset "Core verified spectra against the physical-frame Hamiltonian" begin
    J = 7 // 2
    x = [0.0 0.5; 0.5 0.0]
    a = [0.0 0.5; -0.5 0.0]
    z = [0.5 0; 0 -0.5]
    bond = kron(x, x) - kron(a, a) + kron(z, z)
    for L in (2, 3, 4)
        # Independent dense ED of the physical ladder: site 2j+1 (2j+2) is leg 1 (2) of rung j.
        bonds = [(2j + 1, 2j + 2, Float64(J) * bond) for j in 0:L-1]
        append!(bonds, [(2j + b, 2mod(j + 1, L) + b, bond) for j in 0:L-1 for b in 1:2])
        H = Matrix(sector_hamiltonian(Chain(2, _ -> bonds), 2L, collect(0:4^L-1)))
        E = eigvals(Symmetric(H))
        spec = verified_ring_spectrum(HL.ladder_model(J), L)
        @test length(levels(spec)) == 4^L
        @test maximum(abs, levels(spec) - E) < 1e-12
        @test maximum(s.radius for s in spec) < 1 // 10^10
        # The physical frame (no RP, no flip) has the same spectrum.
        @test maximum(abs, levels(verified_ring_spectrum(HL.physical_model(J), L)) - E) < 1e-12
        # Exact second trace moment; L = 2 has doubled periodic leg bonds.
        expected = big(4)^L * (3 // 16) * (L * J^2 + (L == 2 ? 8 : 2L))
        total = sum(s.multiplicity * sum((interval(c) + interval(-1, 1) *
            to_interval(s.radius))^2 for c in s.centers) for s in spec)
        @test HL.Q(inf(total)) <= expected <= HL.Q(sup(total))
    end
end

@testset "Recurrence rejection is a positive proof" begin
    @test HL.rejection(big(1) // 2, big(0) // 1) !== nothing
    @test HL.rejection(big(0) // 1, big(1) // 2) !== nothing
    @test HL.rejection(big(1) // 100, big(1) // 100) === nothing
    @test_throws ErrorException HL.rejection(big(-1) // 100, big(1) // 100)
    # Exact rational inequalities used in the arbitrary-size argument of METHOD.md.
    @test 1 + 1 + 1 // 2 > 2
    @test 4 / (1 + 9 // 4)^2 == 64 // 169 < 1 // 2
end

@testset "Certified coupling continuation" begin
    seed = Dict("n0" => "4", "beta0" => "2", "eps_s" => "3508321//100000000",
        "eps_t" => "123121//5000000")
    c = HL.continuation(seed, big(1) // 10000)
    @test c !== nothing
    @test HL.q(c["eps_s_upper"]) > HL.q(seed["eps_s"])
    @test HL.q(c["eps_t_upper"]) > HL.q(seed["eps_t"])
    @test HL.q(c["gap_lower_bound"]) > 0
end

@testset "Continuous-temperature negative certificate and coverage" begin
    model = HL.ladder_model(2 // 1)
    specs = Dict(L => verified_ring_spectrum(model, L) for L in (2, 4))
    record = HL.exclude_temperatures(specs, 2)
    @test HL.check_exclusion(specs, record)
    bad = deepcopy(record)
    bad["windows"][1]["beta_from"] = "1//100"
    @test_throws ErrorException HL.check_exclusion(specs, bad)
    bad = deepcopy(record)
    bad["tail"]["eps_s_lower"] = "0//1"
    @test_throws ErrorException HL.check_exclusion(specs, bad)
end

@testset "Stored suite and core certificates (no spectra)" begin
    @test HL.check()
    suite = parse_exact_json(read(HL.SUITE, String))
    bad = deepcopy(suite["metadata"])
    bad["source_sha256"][first(keys(bad["source_sha256"]))] = "0"^64
    @test_throws ErrorException verify_model_metadata(bad, HL.BASE; HL.HASHED...)
    record = suite["parameters"][1]
    J = Rational{Int}(HL.q(record["J"]))
    @test HL.check_coupling(J, record, nothing)
    # A tampered upper defect changes the core interval bound; a tampered witness is caught.
    bad = deepcopy(record)
    row = bad["grid"][findfirst(r -> r["status"] == "certified", bad["grid"])]
    row["eps_t_upper"] = string(HL.q(row["eps_t_upper"]) + 1 // 1000)
    @test_throws ErrorException HL.check_coupling(J, bad, nothing)
    bad = deepcopy(record)
    row = bad["grid"][findfirst(r -> r["status"] == "rejected_by_recurrence", bad["grid"])]
    row["rejection"]["step"] = "99"
    @test_throws ErrorException HL.check_coupling(J, bad, nothing)
end
