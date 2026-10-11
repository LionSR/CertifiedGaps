include(joinpath(@__DIR__, "temperature_lift.jl"))
using Test

const lift_data = SC.read_exact_json(joinpath(@__DIR__, "temperature-lift.json"))
const rat = SC.parse_rational
const asint = x -> Int(SC.integer(x))

@testset "Rational powers with independent integer comparisons" begin
    for x in (frac(0), frac(1, 10^8), frac(2, 3), frac(1), frac(7, 3), frac(1234567)),
        r in (frac(0), frac(1), frac(2), frac(5, 3), frac(67, 40), frac(30, 13), frac(60, 13))

        lo, hi = power_interval(x, r)
        a, b = Int(numerator(r)), Int(denominator(r))
        @test 0 <= lo <= hi && lo^b <= x^a <= hi^b
        @test hi - lo <= frac(1, big(10)^24)
    end
    @test_throws AssertionError power_interval(frac(-1), frac(1, 2))
    @test_throws AssertionError power_interval(frac(1), frac(-1))
    for b in 2:15, n in (big(0), big(1), big(2), big(1234)^b - 1, big(1234)^b, big(1234)^b + 1)
        root = integer_root(n, b)
        @test root^b <= n < (root + 1)^b
    end
end

@testset "Finite distributions check the residual and lifting inequalities" begin
    # The complement contains one J list and three identical N lists.
    # Its finite-list value is computed separately from the endpoint formula.
    for w0 in (frac(1), frac(101, 100), frac(11, 10)),
        tail_j in ([frac(1, 100), frac(1, 50)], [frac(0)]),
        tail_n in ([frac(1, 10), frac(1, 25)], [frac(1, 100)])

        tau = w0 + sum(tail_j)
        D = tau + 3sum(tail_n)
        t = frac(99, 100)
        for s in 1:6
            actual = sum(x^s for x in tail_j) + 3sum(x^s for x in tail_n)
            @test actual <= lifted_residual(D, tau, t, frac(s))
        end
    end
    for z in ([frac(1), frac(1, 5), frac(1, 7)], [frac(6, 5), frac(1, 10)]), r in 1:6
        ratios = z[2:end] ./ first(z)
        R = (sum(z) - 1)^r
        @test sum(x^r for x in ratios) <= R
        purity = sum(x^(2r) for x in z) / sum(x^r for x in z)^2
        @test purity >= inv((1 + R)^2)
    end
end

function check_lift_record(record, inputs)
    d = record["initialization"]
    k, trial = asint(d["moment_length"]), asint(d["trial_length"])
    r, h = rat(d["temperature_multiplier"]), rat(d["h"])
    @assert trial == 120 && iseven(k) && 12 <= k <= trial / 2 && r >= 1
    @assert h == frac(trial, 2k)
    moment, twisted = source1_moment_bounds(k, inputs)
    D, E, tau = (rat(d[key]) for key in ("D", "E", "tau"))
    @assert rat(d["moment_upper"]) == moment && rat(d["twisted_moment_upper"]) == twisted
    @assert D^denominator(r) >= moment^numerator(r)
    @assert E^denominator(r) >= twisted^numerator(r)
    @assert tau == (D + 3E) / 4
    t = rat(d["initial_dominant_lower"])
    exponent = 2h - 1
    @assert 0 < t < tau && 3t > D
    @assert D^denominator(exponent) * t^numerator(exponent) < 1
    for row in d["refinements"]
        @assert rat(row["old_lower"]) == t
        residual = lifted_residual(D, tau, t, 2h)
        @assert residual == rat(row["residual"])
        s = rat(row["new_lower"])
        @assert t < s < tau
        # Rearrange the fractional-power inequality without calling its enclosure.
        @assert s^numerator(2h) < (1 - residual)^denominator(2h)
        t = s
    end
    @assert t == rat(d["dominant_lower"])
    residual = rat(d["spatial_residual_upper"])
    dominant = rat(d["spatial_dominant_power_lower"])
    @assert residual == lifted_residual(D, tau, t, h)
    @assert dominant > 0 && dominant^denominator(h) <= t^numerator(h)
    Q = first(source1_moment_bounds(trial, inputs))
    @assert Q == rat(d["physical_partition_upper"]) > 1
    R = rat(d["physical_excited_ratio_upper"])
    @assert R^denominator(r) >= (Q - 1)^numerator(r)
    p, q = rat(d["p"]), rat(d["q"])
    @assert p == 1 - inv((1 + residual / dominant)^2)
    @assert q == 1 - inv((1 + R)^2)
    chain = record["chain"]
    @assert asint(chain["even_length_threshold"]) == 60
    n, beta = 60, frac(49, 4) * r
    @assert rat(chain["initial_beta"]) == beta
    claim = rat(record["uniform_gap_strict_lower"])
    @assert rat(chain["gap_strict_lower"]) == claim
    grid = rat(chain["rounding_grid"])
    for row in chain["finite_intervals"]
        @assert (rat(row["p"]), rat(row["q"])) == (p, q)
        @assert asint(row["length_from"]) == n && asint(row["length_to"]) == 2n
        @assert rat(row["beta"]) == beta
        A = (1 - p)^2 * (1 - q)
        @assert frac(1, 2) < A < 1 && rat(row["interpolated_purity_lower"]) == A
        @assert first(log_enclosure(A / (1 - A))) / beta == rat(row["gap_lower"]) > claim
        pn, qn = rat(row["next_p"]), rat(row["next_q"])
        @assert (inv(A) - 1)^2 / 2 <= pn < (inv(A) - 1)^2 / 2 + grid
        B = (1 - q)^2 * (1 - pn)
        @assert (inv(B) - 1)^2 / 2 <= qn < (inv(B) - 1)^2 / 2 + grid
        p, q = pn, qn
        n *= 2
        beta *= 2
    end
    tail = chain["tail"]
    @assert asint(tail["length_threshold"]) == n && rat(tail["beta"]) == beta
    @assert (rat(tail["p"]), rat(tail["q"])) == (p, q)
    u, C = rat(tail["u"]), rat(tail["C"])
    @assert u == max(p, q) && 0 < u < frac(1, 6)
    @assert sum((1 - u)^(-j) for j in 1:3)^2 / 2 <= C
    @assert C * u < 1 && C * (1 - 3u) > 3
    @assert first(log_enclosure(inv(C * u))) / beta == rat(tail["gap_lower"]) > claim
    @assert asint(record["uniform_even_length_threshold"]) == 36
    @assert frac(7, 100) < claim < rat(lift_data["bridge"]["gap_lower"])
    return true
end

@testset "Saved certificates and input identity" begin
    @test lift_data["status"] == "exact-scalar-checks-passed"
    @test lift_data["script_sha256"] ==
          bytes2hex(sha256(read(joinpath(@__DIR__, "temperature_lift.jl"))))
    @test lift_data["bootstrap_script_sha256"] ==
          bytes2hex(sha256(read(joinpath(@__DIR__, "optimize_bounds.jl"))))
    @test lift_data["reproduction_verification_sha256"] ==
          bytes2hex(sha256(read(joinpath(@__DIR__, "../results/verification.json"))))
    @test lift_data["source_integrity"] == SC.parse_exact_json(SC.JSON3.write(SC.verify_source()))
    for n in (72, 120)
        path = joinpath(SC.SOURCE, "evidence/trials/certificate$n.json")
        @test lift_data["trial_sha256"][string(n)] == bytes2hex(sha256(read(path)))
    end
    p, q, details = tightened_source1(SC.source1_initialization())
    @test details == lift_data["inputs"]
    bridge = lift_data["bridge"]
    A = (1 - p)^2 * (1 - q)
    @test rat(bridge["p"]) == p && rat(bridge["q"]) == q
    @test rat(bridge["purity_lower"]) == A
    @test rat(bridge["gap_lower"]) == first(log_enclosure(A / (1 - A))) / frac(49, 4)
    @test asint(bridge["length_from"]) == 36 && asint(bridge["length_to"]) >= 60
    inputs = Dict(key => rat(details[key]) for key in ("U", "V", "mJ", "mN"))
    for record in values(lift_data["results"])
        @test check_lift_record(record, inputs)
    end
    for (key, value) in
        (("p", "0"), ("q", "0"), ("dominant_lower", "2"), ("temperature_multiplier", "1"))
        forged = deepcopy(lift_data["results"]["optimized"])
        forged["initialization"][key] = value
        @test_throws AssertionError check_lift_record(forged, inputs)
    end
    forged = deepcopy(lift_data["results"]["optimized"])
    forged["chain"]["finite_intervals"][2]["length_from"] = "122"
    @test_throws AssertionError check_lift_record(forged, inputs)
end
