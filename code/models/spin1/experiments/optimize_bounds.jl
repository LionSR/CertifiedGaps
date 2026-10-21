# Exact scalar improvements using the unchanged, reproduced spin-one certificates.
# Run: julia --project=models/spin1 models/spin1/experiments/optimize_bounds.jl
# No floating-point number is used to accept a mathematical inequality.
import Spin1Certificates as SC
using Dates, SHA
using CertifiedGaps.ExactArithmetic: RationalBig, frac, log_enclosure, first_true
using CertifiedGaps.Recurrence: rational_gap_chain
const GRID = big(10)^512

function squared_filter(coefficients, shift)
    P = fill(frac(0), 2length(coefficients) - 1 + shift)
    for i in eachindex(coefficients), j in eachindex(coefficients)
        P[i+j-1+shift] += coefficients[i] * coefficients[j]
    end
    return P
end
poly(P, x) = sum(c * x^(i - 1) for (i, c) in enumerate(P))

"Check positivity of the full exterior rays via shifted polynomial coefficients."
function check_filter(P, B, cap)
    for sign in (-1, 1)
        shifted = [
            sum(binomial(i, j) * P[i+1] * sign^i * cap^(i - j) for i in j:(length(P)-1)) for
            j in 0:(length(P)-1)
        ]
        @assert first(shifted) > B && all(>=(0), shifted[2:end])
    end
end

function tightened_source1(baseline; readable = false)
    filters = baseline["filters"]
    P1, P2 = (frac.(r["coefficients"]) for r in filters)
    B1, B2 = (SC.parse_rational(r["upper"]) for r in filters)
    if readable
        U, V = frac(1001351, 10^6), frac(871256, 10^6)
    else
        _, U = first_true(x -> poly(P1, x) > B1, frac(1), frac(100139, 100000))
        _, V = first_true(x -> poly(P2, x) > B2, frac(8, 10), frac(872, 1000))
    end
    check_filter(P1, B1, U)
    check_filter(P2, B2, V)
    mJ, mN = frac(2131441, 2000000), frac(556631, 2000000)
    non(k) = 3SC.capped_mass(mN, V^12, div(k, 12))
    if readable
        lo = frac(999997, 10^6)
    else
        lo, _ =
            first_true(x -> SC.capped_mass(mJ, x^12, 6) + non(72) >= 1, frac(999, 1000), frac(1))
    end
    @assert SC.capped_mass(mJ, lo^12, 6) + non(72) < 1
    @assert mJ > lo^12 && 2lo^12 > mJ
    tail(k) = non(k) + (mJ - lo^12)^div(k, 12)
    p = 1 - (1 + tail(36) / lo^36)^(-2)
    q = 1 - (U^72 + tail(72))^(-2)
    return p,
    q,
    Dict(
        "U" => string(U),
        "V" => string(V),
        "dominant_modulus_lower" => string(lo),
        "mJ" => string(mJ),
        "mN" => string(mN),
    )
end

function tightened_source2()
    eps = frac(1, 100000)
    moments = [RationalBig[] for _ in 1:3]
    for (zz, pp, cc) in eachrow(SC.CENTERS2)
        z, p, c = frac.((zz, pp, cc), 10^6)
        push!(moments[1], (z + 3p + 8c) / 12)
        push!(moments[2], (z + 3p - 4c) / 12)
        push!(moments[3], (z - p) / 4)
    end
    caps = RationalBig[]
    for (mom, a, b, lower, upper) in (
        (moments[1], -2, 10, frac(1), frac(100245, 100000)),
        (moments[3], -25, 100, frac(8, 10), frac(88643, 100000)),
    )
        P = squared_filter(frac.([a, 0, b]), 6)
        B = sum(c * m for (c, m) in zip((a^2, 2a*b, b^2), mom)) + eps*(abs(a)+abs(b))^2
        _, cap = first_true(x -> poly(P, x) > B, lower, upper)
        check_filter(P, B, cap)
        push!(caps, cap)
    end
    U, V = caps
    mi, mo, mn = (last(mom) + eps for mom in moments)
    function bounds(k)
        i = SC.capped_mass(mi, U^10, div(k, 10))
        o = mo^div(k, 10)
        n = SC.capped_mass(mn, V^10, div(k, 10))
        return i + 2o + 3n, i + 2o
    end
    m, p = bounds(30)
    q, _ = bounds(120)
    _, D = first_true(x -> x^2 > m^5, frac(1), frac(2))
    _, E = first_true(x -> x^2 > p^5, frac(1), frac(2))
    @assert q > 1
    _, Rq = first_true(x -> x^2 > (q - 1)^5, frac(0), frac(1))
    Ti, T0 = (D + 11E)/12, (D + 3E)/4
    function residual(t)
        @assert 0 < t < Ti < T0 < D
        vertices = ((frac(0), frac(0)), (Ti-t, frac(0)), (Ti-t, T0-Ti), (frac(0), T0-t))
        return maximum(x^2 + y^2/2 + (D-t-x-y)^2/3 for (x, y) in vertices)
    end
    t = frac(88, 100)
    @assert D*t^3 < 1 && 2t > D
    refinements = Any[]
    for _ in 1:20
        R = residual(t)
        # The previous rational lower endpoint ensures strict contradiction.
        s, _ = first_true(x -> x^4 + R^2 >= 1, frac(0), frac(1))
        s <= t && break
        @assert s^4 + R^2 < 1
        push!(
            refinements,
            Dict("old_lower" => string(t), "new_lower" => string(s), "R" => string(R)),
        )
        t = s
    end
    spatial = 1 - (1 + residual(t)/t^2)^(-2)
    temporal = 1 - (1 + Rq)^(-2)
    return spatial,
    temporal,
    Dict(
        "U" => string(U),
        "V" => string(V),
        "mi" => string(mi),
        "mo" => string(mo),
        "mn" => string(mn),
        "D" => string(D),
        "E" => string(E),
        "Ti" => string(Ti),
        "T0" => string(T0),
        "Rq" => string(Rq),
        "dominant_weight_lower" => string(t),
        "refinements" => refinements,
    )
end

function main(output = joinpath(@__DIR__, "bounds.json"))
    SC.write_result(output, Dict("status" => "running"))
    started = time()
    provenance = SC.verify_source()
    reproduction_path = joinpath(@__DIR__, "../results/verification.json")
    reproduction = SC.read_exact_json(reproduction_path)
    @assert reproduction["status"] == "passed"
    @assert reproduction["implementation_sha256"] == SC.implementation_fingerprint()
    for (name, component) in reproduction["components"]
        path = joinpath(@__DIR__, "../results", "$name.json")
        @assert bytes2hex(sha256(read(path))) == component["sha256"]
    end
    baseline1, baseline2 = SC.source1_initialization(), SC.source2_initialization()
    p1, q1 = SC.parse_rational.(baseline1["initial_defects"])
    p2 = 1 - SC.parse_rational(baseline2["spatial_purity_lower"])
    q2 = 1 - SC.parse_rational(baseline2["temporal_purity_lower"])
    t1p, t1q, details1 = tightened_source1(baseline1)
    t2p, t2q, details2 = tightened_source2()
    @assert t1p < p1 && t1q < q1 && t2p < p2 && t2q < q2
    records = Dict{String,Any}()
    for (name, p, q, n, beta, details) in (
        ("original_source1_seed", p1, q1, 36, frac(49, 4), Dict()),
        ("original_source2_seed", p2, q2, 60, frac(105, 4), Dict()),
        ("tightened_source1_seed", t1p, t1q, 36, frac(49, 4), details1),
        ("tightened_source2_seed", t2p, t2q, 60, frac(105, 4), details2),
    )
        record = rational_gap_chain(p, q, n, beta; grid = GRID)
        record["initialization"] = details
        records[name] = record
        println("$name: gamma > $(record["gap_strict_lower_decimal"]) for every even L >= $n")
        flush(stdout)
    end
    rp, rq, rd = tightened_source1(baseline1; readable = true)
    readable =
        rational_gap_chain(rp, rq, 36, frac(49, 4); steps = 9, grid = big(10)^100, tail_C = frac(5))
    readable["initialization"] = rd
    @assert SC.parse_rational(readable["gap_strict_lower"]) > frac(1, 60)
    readable["simple_gap_strict_lower"] = "1//60"
    records["readable_source1_seed"] = readable
    println("readable_source1_seed: gamma > 1/60 for every even L >= 36")
    SC.write_result(
        output,
        Dict(
            "status" => "exact-scalar-checks-passed",
            "results" => records,
            "source_integrity" => provenance,
            "reproduction_verification_sha256" => bytes2hex(sha256(read(reproduction_path))),
            "upstream_commit" => SC.UPSTREAM_COMMIT,
            "execution_mode" => "fresh-scalar-calculation-on-reproduced-finite-inputs",
            "new_thermal_or_trial_contractions" => false,
            "rounding_grid" => string(inv(frac(GRID))),
            "script_sha256" => bytes2hex(sha256(read(@__FILE__))),
            "julia" => string(VERSION),
            "completed_utc" => string(now(UTC)),
            "wall_seconds" => time()-started,
            "scope" => "Exact finite arithmetic plus the cited paper's analytic transfer, interpolation, and bootstrap arguments; not a formal proof of those arguments.",
        ),
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(isempty(ARGS) ? joinpath(@__DIR__, "bounds.json") : only(ARGS))
end
