# Rational-only ports of scalars/check-{moment-initializations,bootstrap-scalars}.py.
"The core's `squaring_defect`, restricted to the update's domain 0 <= u < 1."
function purity_update(u)
    @assert 0 <= u < 1
    return squaring_defect(u)
end
function capped_mass(m, cap, exponent)
    @assert m >= 0 && cap > 0 && exponent >= 1
    count = fld(m, cap)
    return count * cap^exponent + (m - count * cap)^exponent
end
function ceiling_grid(x)
    return frac(cld(numerator(x) * 10^7, denominator(x)), 10^7)
end

function source1_initialization()
    eps = frac(30, 10^6)
    moments = [
        vcat(fill(frac(0), 4), [frac(z + k * p, 4 * 10^6) for (z, p) in eachrow(CENTERS1)]) for
        k in (3, -1)
    ]
    caps = [frac(100139, 100000), frac(872, 1000)]
    quartic = [7, 17, -280, -185, 1000]
    for (sign, expected) in
        ((1, [559, 2902, 5165, 3815, 1000]), (-1, [895, 3978, 6275, 4185, 1000]))
        actual = [sum(binomial(i, j) * quartic[i+1] * sign^i for i in j:4) for j in 0:4]
        @assert actual == expected && minimum(actual) > 0
    end
    @assert caps[1] > 1 && 12 - 394caps[2]^2 + 1000caps[2]^4 > 0 && caps[2]^2 > frac(197, 1000)
    filters = Any[]
    for (r, coefficients) in enumerate((quartic, [12, 0, -394, 0, 1000]))
        P = zeros(Int, 13)
        for i in 0:4, j in 0:4
            P[i+j+5] += coefficients[i+1] * coefficients[j+1]
        end
        upper = sum(c * m + abs(c) * eps for (c, m) in zip(P, moments[r]))
        margins = RationalBig[]
        for sign in (-1, 1)
            shifted = [
                sum(binomial(i, j) * P[i+1] * sign^i * caps[r]^(i - j) for i in j:12) for j in 0:12
            ]
            @assert shifted[1] > upper && all(>=(0), shifted[2:end])
            push!(margins, shifted[1] - upper)
        end
        push!(
            filters,
            Dict("coefficients" => P, "upper" => string(upper), "margins" => string.(margins)),
        )
    end
    up = [m[13] + eps for m in moments]
    cap(n, r, z) = capped_mass(up[r], z^12, div(n, 12))
    non(n) = 3cap(n, 2, caps[2])
    lo = frac(999, 1000)
    @assert cap(72, 1, lo) + non(72) < 1
    @assert up[1] > lo^12 && 2lo^12 > up[1]
    tail(n) = non(n) + (up[1] - lo^12)^div(n, 12)
    p = 1 - (1 + tail(36) / lo^36)^(-2)
    q = 1 - (caps[1]^72 + tail(72))^(-2)
    @assert frac(47915, 10^6) < p < frac(47916, 10^6)
    @assert frac(181520, 10^6) < q < frac(181521, 10^6)
    initial = (p, q)
    expected = [
        (frac(151249, 2500000), frac(433453, 2500000)),
        (frac(343303, 5000000), frac(1632381, 10000000)),
        (frac(44601, 625000), frac(361773, 2500000)),
        (frac(632939, 10000000), frac(1055161, 10000000)),
        (frac(375793, 10000000), frac(445941, 10000000)),
        (frac(84513, 10000000), frac(27493, 5000000)),
    ]
    rounds = Any[]
    for pair in expected
        exact_p = purity_update(1 - (1 - p)^2 * (1 - q))
        p = ceiling_grid(exact_p)
        @assert exact_p <= p < exact_p + frac(1, 10^7) && p < 1
        exact_q = purity_update(1 - (1 - q)^2 * (1 - p))
        q = ceiling_grid(exact_q)
        @assert exact_q <= q < exact_q + frac(1, 10^7) && q < 1
        @assert (p, q) == pair
        push!(rounds, string.([p, q]))
    end
    @assert max(p, q) < frac(1, 100)
    @assert 36 * 2^6 == 2304 && frac(49, 4) * 2^6 == 784
    return Dict(
        "filters" => filters,
        "initial_defects" => string.([initial...]),
        "rounded_defects" => rounds,
        "spatial_length" => 2304,
        "temporal_length" => 4608,
        "beta" => "784",
        "purity_lower" => "99/100",
    )
end

function source2_initialization()
    eps = frac(1, 100000)
    invariant, conjugate, nontrivial = RationalBig[], RationalBig[], RationalBig[]
    for (zz, pp, cc) in eachrow(CENTERS2)
        z, p, c = frac.((zz, pp, cc), 10^6)
        push!(invariant, (z + 3p + 8c) / 12)
        push!(conjugate, (z + 3p - 4c) / 12)
        push!(nontrivial, (z - p) / 4)
    end
    U, V = frac(100245, 100000), frac(88643, 100000)
    mi, mo, mn = frac(1042655, 10^6), frac(40041, 10^6), frac(401625, 10^6)
    for (mom, m) in ((invariant, mi), (conjugate, mo), (nontrivial, mn))
        @assert m > last(mom) + eps
    end
    for (m, cap) in ((mi, U), (mn, V))
        @assert cap^10 < m < 2cap^10
    end
    for (mom, a, b, v) in ((invariant, -2, 10, U), (nontrivial, -25, 100, V))
        upper = sum(c * m for (c, m) in zip((a^2, 2a * b, b^2), mom)) + eps * (abs(a) + abs(b))^2
        @assert b > 0 && a + b * v^2 > 0
        @assert v^6 * (a + b * v^2)^2 > upper
    end
    function bounds(k)
        ib, ob, nb =
            capped_mass(mi, U^10, div(k, 10)), mo^div(k, 10), capped_mass(mn, V^10, div(k, 10))
        @assert ib == U^k + (mi - U^10)^div(k, 10)
        @assert nb == V^k + (mn - V^10)^div(k, 10)
        return ib + 2ob + 3nb, ib + 2ob
    end
    m, p = bounds(30)
    q, _ = bounds(120)
    D, E, Ti, T0, Rq =
        frac(1450, 1000), frac(1202, 1000), frac(1225, 1000), frac(1264, 1000), frac(684, 10000)
    @assert m^5 < D^2 && p^5 < E^2
    @assert (D + 3E) / 4 <= T0 && (D + 11E) / 12 <= Ti
    @assert q >= 1 && (q - 1)^5 < Rq^2
    @assert D * frac(88, 100)^3 < 1 && 2frac(88, 100) > D
    function R(t)
        @assert t < Ti < T0 < D
        vertices = ((frac(0), frac(0)), (Ti - t, frac(0)), (Ti - t, T0 - Ti), (frac(0), T0 - t))
        return maximum(x^2 + y^2 / 2 + (D - t - x - y)^2 / 3 for (x, y) in vertices)
    end
    for (t, s) in ((frac(88, 100), frac(99, 100)), (frac(99, 100), frac(998, 1000)))
        @assert s^4 + R(t)^2 < 1
    end
    S, T = (1 + R(frac(998, 1000)) / frac(998, 1000)^2)^(-2), (1 + Rq)^(-2)
    @assert S > frac(7, 8) && T > frac(7, 8)
    return Dict(
        "spatial_purity_lower" => string(S),
        "temporal_purity_lower" => string(T),
        "spatial_length" => 60,
        "temporal_length" => 120,
        "beta" => "105/4",
        "purity_lower" => "7/8",
    )
end

function bootstrap_checks()
    records = Any[]
    for (name, u0, C, n0, beta) in (
        ("source1", frac(1, 100), frac(5), 2304, frac(784)),
        ("source2", frac(1, 8), frac(79, 10), 60, frac(105, 4)),
    )
        G = bootstrap_coefficient(u0)
        @assert 0 < u0 < frac(1, 6) && G < C && C * u0 < 1
        @assert C * (1 - 3u0) > 3 && (1 - u0)^3 > frac(1, 2)
        for j in 0:6
            u, next_u = (C * u0)^(2^j) / C, (C * u0)^(2^(j + 1)) / C
            @assert next_u == C * u^2 && 0 < next_u <= u <= u0
            pnext = purity_update(1 - (1 - u)^3)
            qnext = purity_update(1 - (1 - u)^2 * (1 - pnext))
            @assert pnext <= next_u && qnext <= next_u
        end
        push!(
            records,
            Dict(
                "source" => name,
                "u0" => string(u0),
                "C" => string(C),
                "coefficient" => string(G),
                "coefficient_margin" => string(C - G),
                "length_threshold" => n0,
                "beta" => string(beta),
                "Cu0" => string(C * u0),
                "ratio_factor" => string(3 / (C * (1 - 3u0))),
            ),
        )
    end
    @assert bootstrap_coefficient(frac(1, 8)) == frac(913952, 117649)
    @assert frac(79, 10) - bootstrap_coefficient(frac(1, 8)) == frac(154751, 1176490)
    @assert frac(9, 2) / frac(97, 100)^2 < 5
    @assert frac(3, 5) / (1 - frac(3, 5) * frac(1, 20)) == frac(60, 97) < 1
    @assert 3 / (frac(79, 10) * (1 - frac(3, 8))) == frac(48, 79) < 1
    return records
end

function run_scalars(output)
    provenance = verify_source()
    result = Dict(
        "status" => "reproduced",
        "source_integrity" => provenance,
        "method" => "exact Rational{BigInt} finite scalar checks",
        "completed_utc" => string(now(UTC)),
        "source1" => source1_initialization(),
        "source2" => source2_initialization(),
        "bootstrap" => bootstrap_checks(),
        "scope" => "Finite arithmetic checks; transfer identities and the all-scale induction remain analytic paper arguments.",
    )
    println(
        "PASS both moment initializations, all six source1 rounds, and both bootstrap constants",
    )
    return write_result(output, result)
end
