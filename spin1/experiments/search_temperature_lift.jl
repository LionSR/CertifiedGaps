# Numerical discovery only. Replay any candidate through temperature_lift.jl.
include(joinpath(@__DIR__, "temperature_lift.jl"))

function numerical_gap(p, q, beta; steps = 10)
    g = BigFloat(1)
    for _ in 1:steps
        # Expanded defects avoid cancellation once the purities approach one.
        a = p * (2 - p) + q * (1 - p)^2
        0 < a < BigFloat(0.5) || return nothing
        g = min(g, (log1p(-a) - log(a)) / beta)
        p = a^2 / (2 * (1 - a)^2)
        b = q * (2 - q) + p * (1 - q)^2
        0 < b < 1 || return nothing
        q = b^2 / (2 * (1 - b)^2)
        beta *= 2
    end
    u = max(p, q)
    C = sum((1 - u)^(-j) for j in 1:3)^2 / 2
    0 < u < 1 / 6 && C * u < 1 && C * (1 - 3u) > 3 || return nothing
    return min(g, -log(C * u) / beta)
end

function search_lift(output = joinpath(@__DIR__, "temperature-lift-search.json"))
    started = time()
    _, _, details = tightened_source1(SC.source1_initialization())
    setprecision(256) do
        U, V, mJ, mN = (BigFloat(SC.parse_rational(details[key])) for key in ("U", "V", "mJ", "mN"))
        J(k) = U^k + (mJ - U^12)^(BigFloat(k) / 12)
        N(k) = V^k + (mN - V^12)^(BigFloat(k) / 12)
        Z(k) = J(k) + 3N(k)
        candidates = Any[]
        for k in 12:2:60, numerator_r in 40:120
            r = BigFloat(numerator_r) / 40
            h = BigFloat(60) / k
            D, E = Z(k)^r, J(k)^r
            tau = (D + 3E) / 4
            t = D^(-inv(2h - 1)) * (1 - big"1e-24")
            3t > D && t < tau || continue
            R(t, s) = max(3 * ((D - t) / 3)^s, (tau - t)^s + 3 * ((D - tau) / 3)^s)
            for _ in 1:40
                R(t, 2h) < 1 || break
                s = (1 - R(t, 2h))^inv(2h) * (1 - big"1e-24")
                s <= t && break
                t = s
            end
            0 < t < tau || continue
            spatial_ratio = R(t, h) / t^h
            physical_ratio = (Z(120) - 1)^r
            defect(x) = x * (2 + x) / (1 + x)^2
            p, q = defect(spatial_ratio), defect(physical_ratio)
            g = numerical_gap(p, q, BigFloat(49) / 4 * r)
            isnothing(g) && continue
            push!(
                candidates,
                Dict(
                    "k" => k,
                    "multiplier" => string(frac(numerator_r, 40)),
                    "gap_estimate" => Float64(g),
                    "p" => Float64(p),
                    "q" => Float64(q),
                ),
            )
        end
        sort!(candidates; by = row -> row["gap_estimate"], rev = true)
        SC.write_result(
            output,
            Dict(
                "status" => "numerical-discovery-only",
                "certified" => false,
                "grid" => "even k=12:2:60; r=40/40:1/40:120/40; trial N=120; ten bootstrap steps",
                "viable_candidates" => length(candidates),
                "top_candidates" => first(candidates, 20),
                "script_sha256" => bytes2hex(sha256(read(@__FILE__))),
                "wall_seconds" => time() - started,
            ),
        )
        println("Best numerical candidate: ", first(candidates))
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    search_lift(isempty(ARGS) ? joinpath(@__DIR__, "temperature-lift-search.json") : only(ARGS))
end
