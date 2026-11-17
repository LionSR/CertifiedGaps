# Fast unit and regression tests for the spin-1 small-seed certificate.
#
#     julia --project=code/julia code/models/spin1_small_seed/test/runtests.jl
include(joinpath(@__DIR__, "..", "src", "Spin1SmallSeed.jl"))
using .Spin1SmallSeed
using Test
using CertifiedGaps.ExactArithmetic

@testset "Spin1SmallSeed" begin
    @testset "reflection positivity" begin
        # rp_certificate verifies -h = Σ w C ⊗ R conj(C) R exactly and throws otherwise.
        rp = Spin1SmallSeed.rp_data()
        @test rp.cert isa Any
        @test string.(rp.weights) == ["1//1", "1//2", "1//2"]
    end

    @testset "exact arithmetic helpers" begin
        @test parse_rational("49//2") == frac(49, 2)
        @test parse_rational("-3/4") == frac(-3, 4)
        @test Spin1SmallSeed.decimal_string(frac(39036761, 10^9), 9) == "0.039036761"
        @test Spin1SmallSeed.decimal_string(frac(1, 2), 9) == "0.500000000"
        @test Spin1SmallSeed.zq(frac(1, 10), frac(1, 10)) == 1 - (frac(9, 10))^3
        @test Spin1SmallSeed.fq(frac(1, 4)) == (frac(1, 16)) / (2 * (frac(3, 4))^2)
    end

    @testset "sector inversion and filters" begin
        table = Spin1SmallSeed.stored_table()
        @test length(table) == 18
        moments = Spin1SmallSeed.sector_moments(table)
        @test length(moments) == 2
        @test all(haskey(m, 12) for m in moments)
        # J_4 = (Z_4 + 3 P_4)/4 exactly from the stored centers.
        z4, p4 = table[(4, 0)], table[(4, 1)]
        @test moments[1][4] == (z4 + 3p4) / 4
        @test moments[2][4] == (z4 - p4) / 4
        # The design's filters re-prove their exterior rays and reproduce the stored caps.
        masses = [m[Spin1SmallSeed.NMAX] + e
            for (m, e) in zip(moments, Spin1SmallSeed.SECTOR_ERRORS)]
        filters = [Spin1SmallSeed.exact_filter(m, e, Spin1SmallSeed.FILTER_SHIFT, q)
            for (m, e, q) in zip(moments, Spin1SmallSeed.SECTOR_ERRORS, Spin1SmallSeed.FILTER_Q)]
        @test parse_rational(filters[1]["cap"]) == frac(10018483311, 10^10)
        @test parse_rational(filters[2]["cap"]) == frac(4361438597, 5 * 10^9)
        # A corrupted chain (rate too small to dominate) must fail the validity predicate.
        bad = Dict("s" => [string(frac(1, 10)), string(frac(1, 2))],
            "t" => [string(frac(1, 10)), string(frac(1, 2))],
            "r" => string(frac(1, 2)), "epsilon" => string(frac(1, 2)))
        @test Spin1SmallSeed.chain_valid(bad) == false
    end

    @testset "frozen input audit" begin
        audit = Spin1SmallSeed.audit_spin1_inputs()
        @test audit["rows_audited"] == 10
        @test length(audit["selected"]) == 10
        @test length(audit["hashes"]) == 2
        trial = Spin1SmallSeed.verify_trial()
        @test trial["length"] == 120
        @test parse(BigInt, trial["positive_margin"]) > 0
    end

    @testset "full exact computation matches the recorded bound" begin
        core = Spin1SmallSeed.compute()
        seed, chain = core["seed"], core["chain"]
        @test parse_rational(seed["p"]) ≈ 0.047914768 atol = 1e-9
        @test parse_rational(seed["q"]) ≈ 0.112578548 atol = 1e-9
        @test parse_rational(chain["gap_strict_lower"]) == frac(39036761, 10^9)
        @test Spin1SmallSeed.chain_valid(chain)
        @test seed["diagonal_threshold_passes"] == true
    end
end
