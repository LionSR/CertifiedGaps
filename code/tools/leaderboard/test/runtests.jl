# Tests of the leaderboard validator on small fixtures (good and bad cases), and of the
# repository's records.yml and news.yml against the Lean records fixture.
#
#     julia --project=tools/leaderboard tools/leaderboard/test/runtests.jl

using Test
using LeaderboardCheck
using LeaderboardCheck: check_decimal!, passed, Report, YAML

const REPO = normpath(joinpath(@__DIR__, "..", "..", "..", ".."))
const FIXTURE = joinpath(@__DIR__, "lean_records.fixture.json")

# A repository with one model whose certificate states the seed (6, 97/4), the bound 1/2 and
# the defects 1/100, next to a floating-point diagnostic that claims more, and a one-step entry
# (a bridge seed) at (4, 3).
const CERT = """{"g": "1/4", "certificate": {"n": 6, "beta": "97/4", "bound_exact_binary": "1/2",
    "epsilon_s": "1/100", "epsilon_t": "1/100"},
    "diagnostic": {"n": 6, "beta": "97/4", "bound": 0.9, "epsilon_s": 0.001, "epsilon_t": 0.001},
    "bridge": {"length_from": 4, "length_to": 8, "beta": "3", "p": "1/50", "q": "1/40"}}"""

const HEAD = """
tiers:
  - {key: lean-unconditional, glyph: 'A', label: 'a', meaning: 'a'}
  - {key: lean-conditional, glyph: 'B', label: 'b', meaning: 'b'}
  - {key: certified, glyph: 'C', label: 'c', meaning: 'c'}
  - {key: target, glyph: 'D', label: 'd', meaning: 'd'}
families:
  - {key: toy, name: Toy, hamiltonian: 'H', slug: toy, cells: sites, gap_source: 'ED'}
"""

record(;
    tier = "lean-conditional",
    bound = """{r: "13480", beta: "97/4", decimal: "0.3921"}""",
    lean = """{name: Leaderboard.Records.Toy.p, r: "13480", beta: "97/4", decimal: "0.3921"}""",
    seed = """{n0: 6, beta0: "97/4"}""",
    extra = "",
    contributors = "[LionSR]",
    params = """{g: "1/4"}""",
    certificate = "code/models/toy/results/certificate.json",
) = """
  - id: toy-p
    family: toy
    parameter: 'g = 1/4'
    params: $params
    lengths: {L0: 6}
    bound: $bound
    seed: $seed
    numerical_gap: "1.19"
    tier: $tier
    $(lean === nothing ? "" : "lean: $lean")
    certificate: $certificate
    date: "2026-10-15"
    pr: 34
    contributors: $contributors
$extra"""

lean_json(; kind = "conditional", L0 = 6, r = "13480", seeds = "") = """
[{"name": "Leaderboard.Records.Toy.p", "kind": "$kind", "L0": $L0, "r": "$r", "beta": "97/4"$seeds},
 {"name": "Leaderboard.Records.Toy.extra", "kind": "conditional", "L0": 6, "r": "2", "beta": "1"}]
"""

function check(rec; json = nothing, head = HEAD)
    mktempdir() do dir
        mkpath(joinpath(dir, "code", "models", "toy", "results"))
        write(joinpath(dir, "code", "models", "toy", "check.jl"), "")
        write(joinpath(dir, "code", "models", "toy", "results", "certificate.json"), CERT)
        write(joinpath(dir, "records.yml"), head * "records:\n" * rec)
        lean = json === nothing ? nothing : joinpath(dir, "lean.json")
        json === nothing || write(lean, json)
        return first(check_records(joinpath(dir, "records.yml"); root = dir, lean_records = lean))
    end
end
fails(report, re) = any(e -> occursin(re, e), report.errors)

@testset "LeaderboardCheck" begin
    @testset "exact numbers" begin
        @test parse_exact("13480") == 13480
        @test parse_exact("97//4") == parse_exact("97/4") == 97 // 4
        @test parse_exact("4134e-1") == parse_exact("413.4") == 2067 // 5
        @test_throws ArgumentError parse_exact(0.5)
        @test_throws ArgumentError parse_exact("1/0")
        @test parse_decimal("0.3921") == (3921 // 10000, 4)
        @test_throws ArgumentError parse_decimal(0.3921)
        lo, hi = bound_enclosure(Dict("r" => "20", "beta" => "784"))
        @test lo <= hi && hi - lo < big(1) // big(10)^20 && Float64(lo) ≈ log(20) / 784
        ok(d) = passed(check_decimal!(Report(), "t", d, lo, hi))
        @test ok("0.00382") && ok("0.0038210")
        @test !ok("0.00383")   # rounded up
        @test !ok("0.0037")    # true, but not the floor at four decimals
    end

    @testset "good records" begin
        @test passed(check(record(); json = lean_json()))
        @test passed(check(record()))   # without the Lean file
        @test passed(check(record(contributors = "['OpenAI (preprint)', LionSR]")))
    end

    @testset "bad records" begin
        @test fails(check(record(tier = "gold")), "tier")
        @test fails(check(record(tier = "unreviewed")), "tier")
        @test fails(
            check(record(bound = """{r: "13480", beta: "97/4", decimal: "0.3922"}""")),
            "exceeds",
        )
        @test fails(
            check(record(bound = """{r: "13480", beta: "97/4", decimal: 0.3921}""")),
            "quoted string",
        )
        @test fails(
            check(record(certificate = "code/models/toy/results/missing.json")),
            "does not exist",
        )
        @test fails(check(record(lean = nothing)), "needs a Lean record")
        @test fails(check(record(contributors = "[]")), "nonempty")
        @test fails(check(record(contributors = "['not/a-handle!']")), "GitHub handle")
        @test fails(check(record(); json = lean_json(L0 = 8)), "L0")
        @test fails(check(record(); json = lean_json(r = "13481")), "r = ")
        @test fails(
            check(record(); json = lean_json(kind = "unconditional")),
            "Lean record is unconditional",
        )
        @test fails(check(record() * record()), "duplicate id")
        target = "  - {id: toy-target, family: toy, parameter: 'x', tier: target, bound: {value: \"1\", decimal: \"1.0\"}}\n"
        @test fails(check(target), "a target has no bound")
        @test fails(check(target), "needs an estimate")
        bad_date = "    history:\n      - {date: \"2026-13-01\", decimal: \"0.39\", lengths: x, tier: certified, pr: 1, method: m}\n"
        @test fails(check(record(extra = bad_date)), "yyyy-mm-dd")
        unordered = """
    history:
      - {date: "2026-10-01", decimal: "0.30", lengths: x, tier: certified, pr: 1, method: m}
      - {date: "2026-10-01", decimal: "0.35", lengths: x, tier: certified, pr: 1, method: m}
"""
        @test fails(check(record(extra = unordered)), "newest first")
        @test passed(check(record(extra = replace(unordered, "0.30" => "0.40"))))
    end

    @testset "the headline tier against the Lean bound" begin
        sharper = """{value: "1/2", decimal: "0.5"}"""   # beats log(13480)/(97/4) ≈ 0.3921
        @test passed(check(record(tier = "certified", bound = sharper)))
        @test fails(check(record(bound = sharper)), "weaker than the headline")
        @test fails(check(record(tier = "certified")), "tier must be lean")
        @test passed(check(record(bound = """{value: "39/100", decimal: "0.39"}""")))
    end

    @testset "rows against their certificates" begin
        # the certificate proves 1/2; the diagnostic's 0.9 is not taken for the certificate
        @test fails(
            check(record(tier = "certified", bound = """{value: "3/5", decimal: "0.6"}""")),
            "less than the row's bound",
        )
        @test fails(check(record(seed = """{n0: 6, beta0: "24"}""")), "states no bound for seed")
        @test fails(check(record(seed = """{n0: 8, beta0: "97/4"}""")), "states no bound for seed")
        @test fails(check(record(params = """{g: "1/2"}""")), "states no bound for seed")
        # a core-format certificate states its parameter point in its params record
        core = Dict{String,Any}(
            "format" => "CertifiedGaps ring gap certificate v1",
            "seed" =>
                Dict("n0" => "6", "beta0" => "97//4", "eps_s" => "1//100", "eps_t" => "1//100"),
            "interval_bound_decimal" => 0.5,
            "params" => Dict("g" => "1//4"),
        )
        c = only(LeaderboardCheck.certificate_candidates(core))
        @test LeaderboardCheck.params_agree(c.scope, Dict("g" => "1/4"))
        @test !LeaderboardCheck.params_agree(c.scope, Dict("g" => "1/2"))
        # a conditional Lean record must assume defects at least the certified ones
        seeds(s, t) = """, "n0": 6, "s": "$s", "t": "$t\""""
        @test passed(check(record(); json = lean_json(seeds = seeds("1/100", "1/50"))))
        @test fails(
            check(record(); json = lean_json(seeds = seeds("1/200", "1/50"))),
            "below the certified defect",
        )
        # a Lean tier may exceed the certificate's bound: the bound is the Lean theorem's
        big = """{r: "1000000", beta: "97/4", decimal: "0.5697"}"""
        bigLean = """{name: Leaderboard.Records.Toy.p, r: "1000000", beta: "97/4", decimal: "0.5697"}"""
        @test passed(
            check(
                record(bound = big, lean = bigLean);
                json = lean_json(r = "1000000", seeds = seeds("1/100", "1/50")),
            ),
        )
        # the second seed of a two-seed record must be a one-step entry of the certificate
        seed2(n, s, t) =
            seeds("1/100", "1/50") * """, "seed2": {"n0": $n, "beta": "3", "s": "$s", "t": "$t"}"""
        @test passed(check(record(); json = lean_json(seeds = seed2(4, "1/50", "1/40"))))
        @test fails(
            check(record(); json = lean_json(seeds = seed2(4, "1/60", "1/40"))),
            "second seed: Lean assumes s",
        )
        @test fails(
            check(record(); json = lean_json(seeds = seed2(5, "1/50", "1/40"))),
            "no one-step entry",
        )
        # a record on a parameter interval needs a Lean record on the same interval
        iv(a, b) = seeds("1/100", "1/50") * """, "interval": ["$a", "$b"]"""
        ivrow = record(params = """{g_min: "1/5", g_max: "1/4"}""")
        @test passed(check(ivrow; json = lean_json(seeds = iv("1/5", "1/4"))))
        @test fails(
            check(ivrow; json = lean_json(seeds = iv("1/5", "1/3"))),
            "is not [g_min, g_max]",
        )
        @test fails(
            check(ivrow; json = lean_json(seeds = seeds("1/100", "1/50"))),
            "the Lean record is a point",
        )
        @test fails(
            check(record(); json = lean_json(seeds = iv("1/5", "1/4"))),
            "the row is a point",
        )
    end

    @testset "challenges and news" begin
        ch(
            link,
        ) = "challenges:\n  - {id: toy-ch, family: toy, title: 'T', need: 'N', link: '$link'}\n"
        @test passed(check(record(); head = HEAD * ch("https://example.org/x")))
        @test fails(check(record(); head = HEAD * ch("code/models/toy/nope.md")), "does not exist")
        @test fails(
            check(record(); head = HEAD * replace(ch("https://example.org"), "toy-ch" => "toy-p")),
            "distinct",
        )
        mktempdir() do dir
            path = joinpath(dir, "news.yml")
            write(
                path,
                """
- {date: "2026-10-15", text: 'Lean theorem', link: 'https://example.org/pull/34', prs: [34]}
- {date: "2026-10-16", text: 'Explainer', link: '/explainer/', records: [toy-p]}
""",
            )
            @test passed(check_news(path; record_ids = Set(["toy-p"])))
            @test !passed(check_news(path; record_ids = Set(String[])))
            write(path, "- {date: \"2026-10-32\", text: 'x', link: 'explainer'}\n")
            r = check_news(path)
            @test fails(r, "yyyy-mm-dd") && fails(r, "link must be")
        end
    end

    @testset "repository data" begin
        path = joinpath(REPO, "home_page", "_data", "records.yml")
        report, summary = check_records(path; root = REPO, lean_records = FIXTURE)
        foreach(println, report.errors)
        @test passed(report)
        @test summary["lean-unconditional"] >= 1 && summary["certified"] >= 40
        ids = Set(String[string(r["id"]) for r in YAML.load_file(path)["records"]])
        news = check_news(joinpath(REPO, "home_page", "_data", "news.yml"); record_ids = ids)
        foreach(println, news.errors)
        @test passed(news)
    end
end
