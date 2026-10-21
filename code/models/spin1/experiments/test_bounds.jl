include(joinpath(@__DIR__, "optimize_bounds.jl"))
using Test

const data = SC.read_exact_json(joinpath(@__DIR__, "bounds.json"))
const rat = SC.parse_rational
const int = x -> Int(SC.integer(x))

@testset "Provenance and finite inputs" begin
    @test data["status"] == "exact-scalar-checks-passed"
    @test data["script_sha256"] == bytes2hex(sha256(read(joinpath(@__DIR__, "optimize_bounds.jl"))))
    @test data["reproduction_verification_sha256"] ==
          bytes2hex(sha256(read(joinpath(@__DIR__, "../results/verification.json"))))
    @test data["source_integrity"] == SC.parse_exact_json(SC.JSON3.write(SC.verify_source()))
    for source in 1:2
        stored = SC.read_exact_json(joinpath(@__DIR__, "../results/source$source.json"))
        for c in stored["certificates"]
            rebuilt =
                SC.thermal_certificate(source, int(c["n"]), int(c["twist"]), SC.integer(c["total"]))
            @test all(c[k] == SC.parse_exact_json(SC.JSON3.write(v)) for (k, v) in rebuilt)
        end
    end
    for n in (72, 120)
        c = SC.read_exact_json(joinpath(SC.SOURCE, "evidence/trials/certificate$n.json"))
        U = SC.integer(c[n == 72 ? "one_bond_energy_numerator" : "U"])
        V = SC.integer(c[n == 72 ? "norm" : "V"])
        @test V > 0 && 500000U + 700741V < 0
    end
end

@testset "Log enclosure independent numerical control" begin
    for x in
        (frac(1), frac(2), frac(3, 2), frac(2)^2000, frac(2)^2000*frac(1999, 1000), frac(10)^300)
        lo, hi = log_enclosure(x)
        setprecision(512) do
            value = log(BigFloat(x))
            @test BigFloat(lo) <= value <= BigFloat(hi)
        end
        @test 0 <= hi - lo <= frac(3, big(10)^24)
    end
    @test_throws AssertionError log_enclosure(frac(1, 2))
end

@testset "Every update, length interval, and uniform tail" begin
    for (_, r) in data["results"]
        claim = rat(r["gap_strict_lower"])
        eps = rat(r["rounding_grid"])
        n, beta = int(r["even_length_threshold"]), rat(r["initial_beta"])
        last_pair = nothing
        for row in r["finite_intervals"]
            p, q = rat(row["p"]), rat(row["q"])
            pn, qn = rat(row["next_p"]), rat(row["next_q"])
            @test isnothing(last_pair) || last_pair == (p, q)
            @test int(row["length_from"]) == n && int(row["length_to"]) == 2n
            @test rat(row["beta"]) == beta
            A = (1-p)^2*(1-q)
            @test frac(1, 2) < A < 1 && A == rat(row["interpolated_purity_lower"])
            # Independently rearrange f(1-A) rather than calling purity_update.
            @test (inv(A)-1)^2/2 <= pn < (inv(A)-1)^2/2 + eps
            B = (1-q)^2*(1-pn)
            @test (inv(B)-1)^2/2 <= qn < (inv(B)-1)^2/2 + eps
            lo, hi = log_enclosure(A/(1-A))
            @test rat(row["gap_lower"]) == lo/beta > claim
            @test rat(row["gap_upper_for_formula"]) == hi/beta
            last_pair = (pn, qn)
            n *= 2
            beta *= 2
        end
        tail = r["tail"]
        @test int(tail["length_threshold"]) == n && rat(tail["beta"]) == beta
        @test last_pair == (rat(tail["p"]), rat(tail["q"]))
        u, C = rat(tail["u"]), rat(tail["C"])
        @test u == max(last_pair...) && 0 < u < frac(1, 6)
        @test sum((1-u)^(-j) for j in 1:3)^2/2 <= C
        @test C*u == rat(tail["r"]) < 1 && C*(1-3u) > 3
        lo, hi = log_enclosure(inv(C*u))
        @test rat(tail["gap_lower"]) == lo/beta > claim
        @test rat(tail["gap_upper_for_formula"]) == hi/beta
    end
end

@testset "Seed inequalities and unsafe cutoff rejection" begin
    baseline1, baseline2 = SC.source1_initialization(), SC.source2_initialization()
    for name in ("original_source1_seed", "original_source2_seed")
        row = first(data["results"][name]["finite_intervals"])
        expected =
            name == "original_source1_seed" ? rat.(baseline1["initial_defects"]) :
            [1-rat(baseline2["spatial_purity_lower"]), 1-rat(baseline2["temporal_purity_lower"])]
        @test [rat(row["p"]), rat(row["q"])] == expected
    end
    for name in ("tightened_source1_seed", "readable_source1_seed")
        r = data["results"][name]
        d = r["initialization"]
        U, V, ell, mJ, mN = (rat(d[k]) for k in ("U", "V", "dominant_modulus_lower", "mJ", "mN"))
        for (filter, cap) in zip(baseline1["filters"], (U, V))
            P, B = frac.(filter["coefficients"]), rat(filter["upper"])
            @test isnothing(check_filter(P, B, cap))
            if name == "tightened_source1_seed"
                @test_throws AssertionError check_filter(P, B, cap-frac(1, 10^12))
            end
        end
        # All instances here have one capped mass, so check the expanded formulas.
        @test V^12 < mN < 2V^12 && ell^12 < mJ < 2ell^12
        non(k) = 3*(V^k+(mN-V^12)^div(k, 12))
        @test ell^72+(mJ-ell^12)^6+non(72) < 1
        R1(k) = non(k)+(mJ-ell^12)^div(k, 12)
        row = first(r["finite_intervals"])
        @test rat(row["p"]) == 1-(1+R1(36)/ell^36)^(-2)
        @test rat(row["q"]) == 1-(U^72+R1(72))^(-2)
    end
    r = data["results"]["tightened_source2_seed"]
    d = r["initialization"]
    U, V, mi, mo, mn, D, E, Ti, T0, Rq =
        (rat(d[k]) for k in ("U", "V", "mi", "mo", "mn", "D", "E", "Ti", "T0", "Rq"))
    moments = [RationalBig[] for _ in 1:3]
    for (z, p, c) in eachrow(SC.CENTERS2)
        push!(moments[1], frac(z+3p+8c, 12*10^6))
        push!(moments[2], frac(z+3p-4c, 12*10^6))
        push!(moments[3], frac(z-p, 4*10^6))
    end
    eps = frac(1, 100000)
    @test (mi, mo, mn) == Tuple(last(m)+eps for m in moments)
    for (mom, a, b, cap) in ((moments[1], -2, 10, U), (moments[3], -25, 100, V))
        P = squared_filter(frac.([a, 0, b]), 6)
        B = a^2*mom[1]+2a*b*mom[2]+b^2*mom[3]+eps*(abs(a)+abs(b))^2
        @test isnothing(check_filter(P, B, cap))
        @test_throws AssertionError check_filter(P, B, cap-frac(1, 10^12))
    end
    @test U^10 < mi < 2U^10 && V^10 < mn < 2V^10
    A(k) = U^k+(mi-U^10)^div(k, 10)
    B(k) = mo^div(k, 10)
    N(k) = V^k+(mn-V^10)^div(k, 10)
    @test (A(30)+2B(30)+3N(30))^5 < D^2
    @test (A(30)+2B(30))^5 < E^2
    @test 1 < A(120)+2B(120)+3N(120)
    @test (A(120)+2B(120)+3N(120)-1)^5 < Rq^2
    @test Ti == (D+11E)/12 && T0 == (D+3E)/4
    t = frac(88, 100)
    @test D*t^3 < 1 && 2t > D
    R2(t) =
        max((D-t)^2/3, (Ti-t)^2+(D-Ti)^2/3, (Ti-t)^2+(T0-Ti)^2/2+(D-T0)^2/3, (T0-t)^2/2+(D-T0)^2/3)
    for row in d["refinements"]
        s = rat(row["new_lower"])
        @test rat(row["old_lower"]) == t && rat(row["R"]) == R2(t)
        @test t < s < Ti && s^4+R2(t)^2 < 1
        t = s
    end
    @test t == rat(d["dominant_weight_lower"])
    row = first(r["finite_intervals"])
    @test rat(row["p"]) == 1-(1+R2(t)/t^2)^(-2)
    @test rat(row["q"]) == 1-(1+Rq)^(-2)
    @test rat(data["results"]["readable_source1_seed"]["gap_strict_lower"]) > frac(1, 60)
end
