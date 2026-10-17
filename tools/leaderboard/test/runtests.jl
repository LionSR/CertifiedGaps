# Fast tests of the leaderboard validator on small fixtures (good and bad cases), plus the real
# records.yml against the Lean records fixture.
#
#     julia --project=tools/leaderboard tools/leaderboard/test/runtests.jl

using Test
using LeaderboardCheck
using LeaderboardCheck: check_decimal!, passed

const REPO = normpath(joinpath(@__DIR__, "..", "..", ".."))
const FIXTURE = joinpath(@__DIR__, "lean_records.fixture.json")

# A minimal repository: one model directory with a certificate and a checker.
function fake_repo(dir)
    mkpath(joinpath(dir, "models", "toy", "results"))
    write(joinpath(dir, "models", "toy", "check.jl"), "")
    write(joinpath(dir, "models", "toy", "results", "certificate.json"), "{}")
    return dir
end

const FAMILY = """
families:
  - key: toy
    name: Toy chain
    hamiltonian: '\$H = -\\sum Z Z\$'
    slug: toy
    cells: sites
"""

record(; tier = "lean-conditional", bound = """{r: "13480", beta: "97/4", decimal: "0.3921"}""",
    lean = """{name: Leaderboard.Records.Toy.p, L0: 6, r: "13480", beta: "97/4", decimal: "0.3921"}""",
    extra = "", certificate = "models/toy/results/certificate.json") = """
  - id: toy-p
    family: toy
    parameter: '\$g = 1/4\$'
    params: {g: "1/4"}
    lengths: {L0: 6, parity: even, unit: sites}
    bound: $bound
    seed: {n0: 6, beta0: "97/4"}
    numerical_gap: {value: "1.19", source: "ED"}
    tier: $tier
    $(lean === nothing ? "" : "lean: $lean")
    certificate: $certificate
    check: models/toy/check.jl
    date: "2026-10-15"
    pr: 34
$extra"""

lean_json(; kind = "conditional", L0 = 6, r = "13480", beta = "97/4",
    hyp = "\"1 - e ≤ P\"", axioms = """["propext", "Classical.choice", "Quot.sound"]""") = """
[{"name": "Leaderboard.Records.Toy.p", "kind": "$kind", "spectrum": "toy.spectrum L",
  "L0": $L0, "r": "$r", "beta": "$beta", "hypothesis": $hyp, "axioms": $axioms,
  "module": "Leaderboard.Records.Toy"},
 {"name": "Leaderboard.Records.Toy.extra", "kind": "conditional", "spectrum": "s", "L0": 6,
  "r": "2", "beta": "1", "hypothesis": "h", "axioms": [], "module": "Leaderboard.Records.Toy"}]
"""

function run_check(yaml; json = nothing)
    mktempdir() do dir
        fake_repo(dir)
        path = joinpath(dir, "records.yml")
        write(path, yaml)
        lean = nothing
        if json !== nothing
            lean = joinpath(dir, "lean_records.json")
            write(lean, json)
        end
        report, summary = check_records(path; root = dir, lean_records = lean)
        return report, summary
    end
end

errors_match(report, re) = any(e -> occursin(re, e), report.errors)

@testset "LeaderboardCheck" begin
    @testset "exact parsing" begin
        @test parse_exact("13480") == 13480
        @test parse_exact("97/4") == 97 // 4
        @test parse_exact("97//4") == 97 // 4
        @test parse_exact("4134e-1") == 2067 // 5
        @test parse_exact("413.4") == 2067 // 5
        @test parse_exact("-3/4") == -3 // 4
        @test parse_exact(6) == 6
        @test_throws ArgumentError parse_exact(0.5)
        @test_throws ArgumentError parse_exact("1/0")
        @test_throws ArgumentError parse_exact("log 2")
        @test parse_decimal("0.3921") == (3921 // 10000, 4)
        @test_throws ArgumentError parse_decimal(0.3921)
        @test_throws ArgumentError parse_decimal("1e-3")
    end

    @testset "decimals rounded down" begin
        lo, hi = bound_enclosure(Dict("r" => "20", "beta" => "784"))
        @test lo <= hi && hi - lo < big(1) // big(10)^20
        @test Float64(lo) ≈ log(20) / 784
        ok(d) = passed(check_decimal!(LeaderboardCheck.Report(), "t", d, lo, hi))
        @test ok("0.00382")
        @test ok("0.0038210")
        @test !ok("0.00383")       # rounded up
        @test !ok("0.0037")        # true but not the floor at four decimals
        v = 1148261 // 1000000
        @test passed(check_decimal!(LeaderboardCheck.Report(), "t", "1.148261", v, v))
        @test !passed(check_decimal!(LeaderboardCheck.Report(), "t", "1.148262", v, v))
    end

    @testset "good fixture" begin
        report, summary = run_check(FAMILY * "records:\n" * record(); json = lean_json())
        @test passed(report)
        @test summary["lean-conditional"] == 1
        @test any(w -> occursin("Toy.extra", w), report.warnings)
        # without the Lean file, the YAML alone validates
        report, _ = run_check(FAMILY * "records:\n" * record())
        @test passed(report)
    end

    @testset "bad fixtures" begin
        bad(yaml; json = nothing) = first(run_check(FAMILY * "records:\n" * yaml; json))
        @test errors_match(bad(record(tier = "gold")), "tier")
        @test errors_match(bad(record(bound = """{r: "13480", beta: "97/4", decimal: "0.3922"}""")),
            "exceeds")
        @test errors_match(bad(record(bound = """{r: "13480", beta: "97/4", decimal: 0.3921}""")),
            "quoted string")
        @test errors_match(bad(record(bound = """{value: "1/2", r: "2", beta: "1", decimal: "0.5"}""")),
            "not both")
        @test errors_match(bad(record(certificate = "models/toy/results/missing.json")),
            "does not exist")
        @test errors_match(bad(record(certificate = "../outside.json")), "relative")
        @test errors_match(bad(record(lean = nothing)), "needs a lean record")
        @test errors_match(bad(record(tier = "certified")), "tier must be lean")
        @test errors_match(bad(record(); json = lean_json(kind = "unconditional", hyp = "null")),
            "Lean record is unconditional")
        @test errors_match(bad(record(); json = lean_json(L0 = 8)), "L0")
        @test errors_match(bad(record(); json = lean_json(r = "13481")), "r = ")
        @test errors_match(bad(record(); json = lean_json(beta = "24")), "beta = ")
        @test errors_match(bad(record(); json = lean_json(axioms = """["sorryAx"]""")),
            "nonstandard axioms")
        @test errors_match(bad(record(); json = lean_json(hyp = "null")), "no hypothesis")
        @test errors_match(bad(record(lean = """{name: Other.p, L0: 6, r: "13480", beta: "97/4"}""")),
            "Leaderboard.Records")
        @test errors_match(bad(record(lean = """{name: Leaderboard.Records.Toy.q, L0: 6, r: "13480", beta: "97/4"}""");
            json = lean_json()), "not in lean_records.json")
        @test errors_match(bad(record() * replace(record(), "toy-p" => "toy-p")), "duplicate id")
        @test errors_match(bad(replace(record(), "family: toy" => "family: nope")), "unknown family")
        target = """
  - id: toy-target
    family: toy
    parameter: 'x'
    bound: {value: "1", decimal: "1.0"}
    tier: target
"""
        @test errors_match(bad(target), "a target has no bound")
        @test errors_match(bad(target), "needs an estimate")
        history = """
    history:
      - {date: "2026-13-01", decimal: "0.39", method: "m", tier: certified}
"""
        @test errors_match(bad(record(extra = history)), "yyyy-mm-dd")
        unordered = """
    history:
      - {date: "2026-10-01", decimal: "0.30", method: "m", tier: certified}
      - {date: "2026-10-01", decimal: "0.35", method: "m", tier: certified}
"""
        @test errors_match(bad(record(extra = unordered)), "newest first")
        older_later = replace(unordered, "\"2026-10-01\", decimal: \"0.35" => "\"2026-10-02\", decimal: \"0.35")
        @test errors_match(bad(record(extra = older_later)), "newest first")
        ordered = replace(unordered, "0.30" => "0.40")
        @test passed(first(run_check(FAMILY * "records:\n" * record(extra = ordered))))
    end

    @testset "headline tier against the Lean bound" begin
        check(yaml; json = lean_json()) = first(run_check(FAMILY * "records:\n" * yaml; json))
        sharper = """{value: "1/2", decimal: "0.5"}"""   # beats log(13480)/(97/4) ≈ 0.3921
        # a certificate sharper than the Lean bound: the headline tier is the certificate's
        @test passed(check(record(tier = "certified", bound = sharper)))
        @test errors_match(check(record(bound = sharper)), "weaker than the headline")
        # a Lean record with the same bound as the headline: the tier must be Lean's
        @test errors_match(check(record(tier = "certified")), "tier must be lean")
        # a smaller headline value is implied by the Lean theorem
        @test passed(check(record(bound = """{value: "39/100", decimal: "0.39"}""")))
        # the Lean theorem must also cover the headline lengths
        fewer = """{name: Leaderboard.Records.Toy.p, L0: 8, r: "13480", beta: "97/4"}"""
        @test errors_match(check(record(lean = fewer)), "weaker than the headline")
        # the kind must match when the tier is Lean's
        @test errors_match(check(record(tier = "lean-unconditional")), "Lean record is conditional")
        # long exact values need a longer rounded-down decimal
        long = """{value: "7066853945252893/18014398509481984", decimal: "0.3922"}"""
        @test errors_match(check(record(tier = "certified", bound = long)), "add bound.precise")
        withp = replace(long, "}" => """, precise: "0.392289198084"}""")
        @test passed(check(record(tier = "certified", bound = withp)))
        wrongp = replace(withp, "0.392289198084" => "0.392289198085")
        @test !passed(check(record(tier = "certified", bound = wrongp)))
    end

    @testset "Lean source scan" begin
        mktempdir() do dir
            lb = joinpath(dir, "leaderboard")
            mkpath(joinpath(lb, "Records"))
            write(joinpath(lb, "Clean.lean"), """
                /-! The audit rejects `sorry` and `native_decide`; this docstring mentions them. -/
                -- an axiom in a comment is fine: axiom foo : False
                theorem t : 1 = 1 := rfl
                def s := "sorry inside a string is fine"
                """)
            @test isempty(scan_lean_sources(lb))
            write(joinpath(lb, "Records", "Bad.lean"), """
                theorem a : 1 = 2 := by sorry
                theorem b : 2 + 2 = 4 := by native_decide
                axiom c : False
                set_option debug.skipKernelTC true in
                @[implemented_by foo] def d := 1
                /- nested /- comment -/ still a comment: sorry -/
                theorem e : True := trivial
                """)
            hits = scan_lean_sources(lb)
            @test length(hits) == 5
            @test any(h -> occursin("Bad.lean:1: sorry", h), hits)
            @test any(h -> occursin("Bad.lean:3: axiom", h), hits)
            @test !any(h -> occursin(":6:", h), hits)
        end
        @test isempty(scan_lean_sources(joinpath(tempdir(), "no-such-leaderboard-dir")))
    end

    @testset "repository records.yml" begin
        report, summary = check_records(joinpath(REPO, "home_page", "_data", "records.yml");
            root = REPO, lean_records = FIXTURE)
        foreach(println, report.errors)
        @test passed(report)
        @test isempty(report.warnings)
        @test summary["lean-unconditional"] >= 1
        @test sum(values(summary)) >= 60
    end
end
