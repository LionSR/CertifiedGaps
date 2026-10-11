using Test, Random
import Spin1Certificates as SC

@testset "exact certificate parsing and signed floor" begin
    mktemp() do path, io
        write(
            io,
            """{"large":1234567890123456789012345678901234567890,"negative":-98765432101234567890,"label":"quote \\\" 123 and 4e2","float":0.5,"flag":true}""",
        )
        close(io)
        value = SC.read_exact_json(path)
        @test SC.integer(value["large"]) ==
              parse(BigInt, "1234567890123456789012345678901234567890")
        @test SC.integer(value["negative"]) == -parse(BigInt, "98765432101234567890")
        @test value["label"] == "quote \" 123 and 4e2"
        @test value["float"] == 0.5 && value["flag"]
        @test_throws ErrorException SC.integer(value["float"])
    end
    @test SC.parse_rational("-3/4") == SC.parse_rational("-3//4") == -big(3)//4
    @test SC.parse_rational(Dict("numerator" => "12", "denominator" => "16")) == big(3)//4
    @test all((x >> k) == fld(x, 2^k) for x in -100:100 for k in (4, 5))
    @test (-1 >> 5) == -1 && div(-1, 32) == 0
end

@testset "source integrity and independent physics" begin
    @test SC.verify_source()["verified_files"] == 43
    @test length(SC.matrix_and_symmetry_checks()) == 12
    @test length(SC.triad_checks()) == 27
    @test length(SC.error_bound_checks()) == 9
    labels, A = SC.trial_matrices()
    @test length(labels) == 26
    @test SC.physical_trial_check(A, SC.trial_blocks(labels, A))["physical_and_transfer_agree"]
end

@testset "fresh thermal recurrence and batch invariance" begin
    for (source, n, twist) in ((1, 4, 0), (1, 4, 1), (2, 6, 2))
        results = [SC.thermal_sector(source, n, twist, m; batch = 7) for m in 0:n]
        total = sum(x.total for x in results)
        @test sum(x.multiplicity for x in results) == 3^n
        @test SC.thermal_certificate(source, n, twist, total)["matches_upstream_exactly"]
        # A different lane width must not change even one integer.
        @test results[1].total == SC.thermal_sector(source, n, twist, 0; batch = 1).total
        @test_throws AssertionError SC.thermal_certificate(source, n, twist, total + 1)
    end
end

@testset "exact moment filters and bootstrap" begin
    @test SC.source1_initialization()["purity_lower"] == "99/100"
    @test SC.source2_initialization()["purity_lower"] == "7/8"
    @test length(SC.bootstrap_checks()) == 2
    @test SC.ceiling_grid(SC.frac(1, 3)) == SC.frac(3333334, 10^7)
    @test_throws AssertionError SC.purity_update(SC.frac(1))
end

@testset "optimized kernels preserve every signed integer coordinate" begin
    rng = MersenneTwister(1234)
    basis = SC.magnetization_basis(Val(4), 0)
    u = rand(rng, (-(Int64(1)<<55)):(Int64(1)<<55), 7, length(basis))
    v = rand(rng, (-(Int64(1)<<55)):(Int64(1)<<55), size(u))
    for source in (1, 2), twist in (source == 1 ? (0:1) : (0:2))
        rows = SC.thermal_rows(basis, source, twist)
        reference, actual = similar(u), similar(u)
        shift = source == 1 ? 5 : 4
        if source == 2 && twist == 2
            refv, gotv = similar(u), similar(u)
            SC.multiply_eisenstein_reference!(reference, refv, rows, u, v)
            SC.multiply_eisenstein!(actual, gotv, rows, u, v, shift)
            @test gotv == (refv .>> shift)
        else
            SC.multiply_real_reference!(reference, rows, u)
            SC.multiply_real!(actual, rows, u, shift)
        end
        @test actual == (reference .>> shift)
    end
end

@testset "exact SU(2) reduction and native block construction" begin
    labels, A = SC.trial_matrices()
    reference = SC.trial_blocks(labels, A)
    native = SC.trial_blocks(labels, A; native = true)
    @test all(x.T == y.T && x.O == y.O for (x, y) in zip(reference, native))
    reduced = SC.su2_trial_blocks(labels, A, native)
    @test [b.dimension for b in reduced] == [26, 56, 48, 24, 7, 1]
    @test sum((2b.spin + 1) * b.dimension for b in reduced) == 676
    for n in (4, 6)
        @test SC.contract_trial_su2(reduced, n) == SC.contract_trial(reference, n)
    end
end
