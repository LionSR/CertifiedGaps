# Independent arithmetic audit of the J⊥ = 4 seed in results/certificate.json (about 40 s).
#
# The spectra are recomputed (no stored radius is used), then the partition sums are redone in
# 256-bit MPFR with directed rounding instead of IntervalArithmetic, and the first two exact
# trace moments Tr H = 0 and Tr H² are checked against the eigenvalue enclosures. Writes
# results/audit.json.
#
#     OPENBLAS_NUM_THREADS=1 julia --project=models/heisenberg_ladder \
#         models/heisenberg_ladder/audit.jl
include(joinpath(@__DIR__, "src", "HeisenbergLadder.jl"))
using .HeisenbergLadder: Q, q, CERTIFICATE, exact_rp, provenance, verified_spectra
using CertifiedGaps: certified_gap_bound
using CertifiedGaps.Results: write_result, parse_exact_json
using Test

"Directed-rounding MPFR enclosure `(lo, hi)` of log Z(β) from every eigenvalue enclosure."
function mpfr_logZ(s, beta)
    setprecision(BigFloat, 256) do
        B = BigFloat(beta); @assert Q(B) == beta
        lo = BigFloat(0); hi = BigFloat(0)
        for (e, r) in zip(s.energies, s.radii)
            lower = setrounding(() -> BigFloat(e) - BigFloat(r), BigFloat, RoundDown)
            upper = setrounding(() -> BigFloat(e) + BigFloat(r), BigFloat, RoundUp)
            lo = setrounding(() -> lo + exp(-B * upper), BigFloat, RoundDown)
            hi = setrounding(() -> hi + exp(-B * lower), BigFloat, RoundUp)
        end
        (setrounding(() -> log(lo), BigFloat, RoundDown),
         setrounding(() -> log(hi), BigFloat, RoundUp))
    end
end

"""
Exact trace-moment checks: Tr H = 0 and Tr H² = 4^L (3/16) L (J² + 2) for L ≥ 3 rungs (three
bonds per rung, each S·S with Tr (S·S)² = 3/16 per state), enclosed by the eigenvalue balls.
"""
function trace_checks(s, L, J)
    low = Q(0); high = Q(0); low2 = Q(0); high2 = Q(0)
    for (e, r) in zip(s.energies, s.radii)
        a = Q(e) - Q(r); b = Q(e) + Q(r)
        low += a; high += b
        low2 += a <= 0 <= b ? 0 : min(a^2, b^2)
        high2 += max(a^2, b^2)
    end
    expected = big(4)^L * (3 // 16) * L * (J^2 + 2)
    @test low <= 0 <= high
    @test low2 <= expected <= high2
    Dict("L" => L, "dimension" => s.dimension, "trace_lower" => string(low),
        "trace_upper" => string(high), "trace2_lower" => string(low2),
        "trace2_upper" => string(high2), "trace2_exact" => string(expected))
end

function audit()
    start = time()
    cert = parse_exact_json(read(CERTIFICATE, String))
    result = only(filter(r -> q(r["J"]) == 4, cert["parameters"]))
    seed = result["best"]
    n = Int(q(seed["n0"])); beta = q(seed["beta0"]); J = 4 // 1
    eps_s, eps_t = q(seed["eps_s_upper"]), q(seed["eps_t_upper"])
    @test n == 4 && beta == 2
    @test exact_rp(J)
    specs = verified_spectra(J)
    a, b = specs[n], specs[2n]
    traces = [trace_checks(a, n, J), trace_checks(b, 2n, J)]
    la = mpfr_logZ(a, beta); lb = mpfr_logZ(b, beta); lc = mpfr_logZ(b, 2beta)
    ps = setrounding(() -> exp(lb[1] - 2la[2]), BigFloat, RoundDown)
    pt = setrounding(() -> exp(lc[1] - 2lb[2]), BigFloat, RoundDown)
    es = setrounding(() -> 1 - ps, BigFloat, RoundUp)
    et = setrounding(() -> 1 - pt, BigFloat, RoundUp)
    @test Q(es) <= eps_s
    @test Q(et) <= eps_t
    @test q(seed["gap_lower_bound"]) <= certified_gap_bound(eps_s, eps_t, beta).bound
    enc(x) = Dict("lo" => string(Q(x[1])), "hi" => string(Q(x[2])))
    write_result(joinpath(@__DIR__, "results", "audit.json"), Dict(
        "metadata" => provenance(Dict("arithmetic" =>
            "256-bit directed MPFR and exact rational trace moments"); certified = true),
        "seed" => Dict("J" => "4//1", "n0" => n, "beta0" => string(beta),
            "eps_s_upper" => string(eps_s), "eps_t_upper" => string(eps_t)),
        "stored_radii_used" => false, "fresh_dimensions" => [a.dimension, b.dimension],
        "traces" => traces, "logZ_n_beta" => enc(la), "logZ_2n_beta" => enc(lb),
        "logZ_2n_2beta" => enc(lc), "eps_s_mpfr_upper" => string(Q(es)),
        "eps_t_mpfr_upper" => string(Q(et)), "seconds" => time() - start,
        "limits" => "Fresh eigensolves share the symmetry builder and norm-bound " *
            "implementation; MPFR partition arithmetic and exact trace moment checks are " *
            "separate. The repository theorem and IEEE/OpenBLAS arithmetic model remain " *
            "assumptions."))
    println("AUDIT PASSED: fresh spectra, exact trace moments, independent directed MPFR ",
        "purities; ", round(time() - start; digits = 1), " s")
end

audit()
