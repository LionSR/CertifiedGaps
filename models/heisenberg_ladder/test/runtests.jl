# Fast tests (under a minute): exact RP, the symmetry builder against an independent full
# Hamiltonian, verified versus numerical spectra, the rejection witness, continuation, the
# all-β exclusion and its checker, and the stored certificate without recomputing spectra.
using Test, LinearAlgebra, IntervalArithmetic
include(joinpath(@__DIR__, "..", "src", "HeisenbergLadder.jl"))
const HL = HeisenbergLadder
using CertifiedGaps.ExactDiagonalization: Chain, sector_hamiltonian
using CertifiedGaps.RingModels: verify_rp, encode_model, decode_model, encode_rp, decode_rp
using CertifiedGaps.Results: parse_exact_json, sha256_file
using CertifiedGaps: certified_gap_bound

@testset "Exact RP and original-to-rotated identities" begin
    for J in HL.COUPLINGS
        @test HL.exact_rp(J)
        model, rp = HL.rp_model(J)
        @test verify_rp(model, rp)
        @test decode_model(encode_model(model)).bond == model.bond
        bad = deepcopy(encode_rp(rp))
        bad["factors"][1]["weight"] = "2//1"
        @test !verify_rp(model, decode_rp(bad))
    end
end

@testset "Sector reconstruction and complete small spectra" begin
    for L in (2, 3, 4)
        J = 7 // 2
        H = zeros(4^L, 4^L)
        for s in 0:4^L-1, (t, h) in HL.transitions(s, L, J)
            H[t+1, s+1] = Float64(h)
        end
        # Independent shared-core graph builder checks the bit transition convention.
        x = [0.0 0.5; 0.5 0.0]
        a = [0.0 0.5; -0.5 0.0]
        z = [0.5 0; 0 -0.5]
        bond = kron(x, x) - kron(a, a) + kron(z, z)
        bonds = [(2j + 1, 2j + 2, Float64(J) * bond) for j in 0:L-1]
        append!(bonds, [(2j + b, 2mod(j + 1, L) + b, bond) for j in 0:L-1 for b in 1:2])
        chain = Chain(2, _ -> bonds)
        @test H == Matrix(sector_hamiltonian(chain, 2L, collect(0:4^L-1)))
        full, rad, _ = HL.verified_eigen(H, 0.0)
        reduced = HL.full_spectrum(L, J; verified = true)
        @test all(abs.(full - reduced.energies) .<= rad .+ reduced.radii)
        # The floating-point scan path must agree with the verified one.
        numerical = HL.full_spectrum(L, J)
        @test all(abs.(numerical.energies - reduced.energies) .<= reduced.radii .+ 1e-10)
        @test sum(d * m for (d, m) in values(HL.sector_dimensions(L))) == 4^L
        # L=2 has doubled periodic leg bonds, so its second trace moment differs.
        expected = big(4)^L * (3 // 16) * (L * J^2 + (L == 2 ? 8 : 2L))
        total = sum((interval(e) + interval(-r, r))^2
            for (e, r) in zip(reduced.energies, reduced.radii))
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
    seed = Dict("n0" => "4//1", "beta0" => "2//1",
        "eps_s_upper" => "3508321//100000000", "eps_t_upper" => "123121//5000000")
    c = HL.continuation(seed, big(1) // 10000)
    @test c !== nothing
    @test HL.q(c["eps_s_upper"]) > HL.q(seed["eps_s_upper"])
    @test HL.q(c["eps_t_upper"]) > HL.q(seed["eps_t_upper"])
    @test HL.q(c["gap_lower_bound"]) > 0
end

@testset "Continuous-temperature negative certificate and coverage" begin
    specs = Dict(L => HL.full_spectrum(L, 2 // 1; verified = true) for L in (2, 4))
    record = HL.exclude_temperatures(specs, 2)
    @test HL.check_exclusion(specs, record)
    bad = deepcopy(record)
    bad["windows"][1]["beta_from"] = "1//100"
    @test_throws ErrorException HL.check_exclusion(specs, bad)
    bad = deepcopy(record)
    bad["tail"]["eps_s_lower"] = "0//1"
    @test_throws ErrorException HL.check_exclusion(specs, bad)
end

@testset "Stored certificate: hashes, recurrence and claims (no spectra)" begin
    cert = parse_exact_json(read(HL.CERTIFICATE, String))
    @test cert["format"] == HL.FORMAT
    # Only model-local sources are enforced; the core hashes are provenance.
    @test Set(String.(keys(cert["metadata"]["source_sha256"]))) == Set(keys(HL.local_hashes()))
    for (file, hash) in cert["metadata"]["source_sha256"]
        @test sha256_file(joinpath(HL.ROOT, file)) == hash
    end
    @test !isempty(cert["metadata"]["core_sha256"])
    bad = deepcopy(cert)
    bad["metadata"]["source_sha256"]["julia/src/Results.jl"] = "0"^64
    @test_throws ErrorException HL.check_header(bad)
    for result in cert["parameters"]
        J = Rational{Int}(HL.q(result["J"]))
        best = result["best"]
        @test best == HL.best_row(result["grid"])
        if best !== nothing
            b = certified_gap_bound(HL.q(best["eps_s_upper"]), HL.q(best["eps_t_upper"]),
                HL.q(best["beta0"])).bound
            @test 0 < HL.q(best["gap_lower_bound"]) <= b
            @test iseven(Int(HL.q(best["n0"])))
        end
        @test HL.check_claim(J, result)
    end
end
