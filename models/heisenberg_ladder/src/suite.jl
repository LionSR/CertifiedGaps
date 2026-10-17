# The certificate suite on top of the core pipeline. Per coupling: verified spectra of the
# rings of 2, 4 and 8 rungs (`verified_ring_spectrum`), every grid seed classified from the
# core seed defects (`seed_bounds`), a core gap certificate (`certify_gap`) at the best
# certified seed, and the model-specific parts: rejection witnesses, all-β exclusions, the
# coupling continuation and the README claims. No numerical scan or older result is an input.

"Seed temperatures β₀ = k/4, k = 1, …, 64."
const BETAS = [big(k) // 4 for k in 1:64]
"Even seed lengths reachable with complete spectra of rings of at most 8 rungs."
const SEED_LENGTHS = (2, 4)
"Rings whose verified spectra the suite needs."
const RINGS = (2, 4, 8)
"""
Couplings in the suite: the five certified ones (four round values and the upper end of the
finite-grid threshold bracket), the two all-β exclusions, and the lower bracket end (proposed
by `scripts/threshold_bisection.jl`).
"""
const COUPLINGS = [4 // 1, 7 // 2, 3 // 1, 5 // 2, 2 // 1, 3 // 2,
    326133 // 131072, 163069 // 65536]
"Couplings at which seeds n₀ ∈ {2, 4} are excluded for every β > 0."
const ALL_BETA_EXCLUDED = (2 // 1, 3 // 2)
"""
README headline claims checked by `check`: a lower bound on the gap from the core
certificate (its interval bound `interval_bound_decimal`), or `nothing` where every grid seed
must be rigorously rejected.
"""
const CLAIMED_GAPS = Dict{Rational{Int},Union{Q,Nothing}}(
    4 // 1 => big(1102404331) // 10^9, 7 // 2 => big(778648930) // 10^9,
    3 // 1 => big(454189766) // 10^9, 5 // 2 => big(42434525) // 10^9,
    2 // 1 => nothing, 3 // 2 => nothing, 326133 // 131072 => nothing,
    163069 // 65536 => big(1540040) // 10^9)
"Coupling half-width of the continuation intervals around certified couplings."
const CONTINUATION_DELTA = big(1) // 10000
"Grid of the lower defect bounds and continuation defects (the core seeds also use 10⁻¹⁵)."
const GRID = big(10)^15
const FORMAT = "heisenberg-ladder suite v2 (core certificates)"
"Model-specific suite record: grid decisions, exclusions, continuations and summaries."
const SUITE = joinpath(BASE, "results", "suite.json")

q(x::AbstractString) = parse(Q, x)
q(x::Integer) = Q(x)
token(x) = string(Q(x))
"File tag of a coupling, e.g. 7//2 -> \"7_2\"."
coupling_tag(J) = denominator(J) == 1 ? string(numerator(J)) :
                  string(numerator(J), "_", denominator(J))
"Core certificate of a certified coupling."
certificate_path(J) = joinpath(BASE, "results", "certificate_J$(coupling_tag(J)).json")

"""
Model-local sources hashed into every result by `Results.model_metadata` and enforced, key set
included, by `Results.verify_model_metadata`: the library and the entry points. The environment
is the core's (`--project=julia`), recorded with the core as provenance only. Scripts,
tests and documentation are not hashed.
"""
const HASHED = (; dirs = ("src",), files = ("certify.jl", "check.jl"))

"Fresh verified spectra (core `BlockSpectrum`s) of the rings of 2, 4 and 8 rungs."
verified_spectra(model) = Dict(L => verified_ring_spectrum(model, L) for L in RINGS)

"Core seed defects of the seed (n, β), reusing the precomputed spectra of rings n and 2n."
seed(model, specs, n, beta) =
    seed_bounds(model, n, beta; spectra = Dict(n => specs[n], 2n => specs[2n]))

"""
    rejection(s, t; steps = 64) -> witness or nothing

Prove from lower defect bounds `(s, t)` that the doubling recurrence cannot converge to zero:
iterate lower bounds (rounded down, which the order-preserving map tolerates) until an
argument z reaches 1/2, after which both components stay ≥ 1/2. `nothing` means no proof,
never a certificate.
"""
function rejection(s::Q, t::Q; steps = 64)
    0 <= s <= 1 && 0 <= t <= 1 || error("invalid lower defects")
    grid = big(10)^24
    for j in 0:steps
        z = 1 - (1 - t) * (1 - s)^2
        # z >= 1/2 implies both next components >= 1/2, hence no attraction to zero.
        z >= 1 // 2 && return Dict("step" => token(j), "z_lower" => token(z))
        sp = round_down(z^2 / (2 * (1 - z)^2), grid)
        zt = 1 - (1 - sp) * (1 - t)^2
        zt >= 1 // 2 && return Dict("step" => token(j), "thermal_z_lower" => token(zt))
        s, t = sp, round_down(zt^2 / (2 * (1 - zt)^2), grid)
    end
    return nothing
end

"Lower end of a defect enclosure, rounded down to `GRID` and clipped at zero."
lower(x::Interval) = max(big(0) // 1, round_down(Q(inf(x)), GRID))

"""
    classify(seed) -> Dict

Classify one grid seed from its core `seed_bounds` result: `certified` (the upper defects
close the core recurrence), `rejected_by_recurrence` (the lower defects carry a `rejection`
witness) or `unresolved`. A failed upper-bound attempt is never called a rejection.
"""
function classify(n, beta, b)
    su, tu, sl, tl = b.eps_s, b.eps_t, lower(b.defect_s), lower(b.defect_t)
    answer = su < 1 && tu < 1 ? certified_gap_bound(su, tu, beta).bound : nothing
    witness = rejection(sl, tl)
    answer !== nothing && witness !== nothing && error("inconsistent classification")
    status = answer !== nothing ? "certified" :
             witness !== nothing ? "rejected_by_recurrence" : "unresolved"
    record = Dict{String,Any}("n0" => token(n), "beta0" => token(beta), "status" => status,
        "eps_s_upper" => token(su), "eps_t_upper" => token(tu),
        "eps_s_lower" => token(sl), "eps_t_lower" => token(tl))
    answer !== nothing && (record["interval_bound"] = token(answer))
    witness !== nothing && (record["rejection"] = witness)
    return record
end

"Highest certified grid row by its interval gap bound, or `nothing`."
function best_row(rows)
    good = filter(r -> r["status"] == "certified", rows)
    return isempty(good) ? nothing : argmax(r -> q(r["interval_bound"]), good)
end

"""
    continuation(seed, delta) -> Dict or nothing

Enlarge a certified core seed `(n0, beta0, eps_s, eps_t)` to every coupling within `delta`.
The centred rung sum satisfies ‖W_L + L/4‖ ≤ L/2, so each log Z_L moves by at most βδL/2 after
a scalar shift that cancels in both purities; this gives s_δ = 1-(1-s)e^{-2βnδ} and
t_δ = 1-(1-t)e^{-4βnδ}, rounded up to `GRID`, and the core interval gap bound.
"""
function continuation(seed, delta)
    n, beta = Int(q(seed["n0"])), q(seed["beta0"])
    s, t = q(seed["eps_s"]), q(seed["eps_t"])
    su = round_up(Q(sup(1 - (1 - to_interval(s)) * exp(-to_interval(2beta * n * delta)))), GRID)
    tu = round_up(Q(sup(1 - (1 - to_interval(t)) * exp(-to_interval(4beta * n * delta)))), GRID)
    b = certified_gap_bound(su, tu, beta).bound
    b === nothing && return nothing
    return Dict("delta" => token(delta), "eps_s_upper" => token(su), "eps_t_upper" => token(tu),
        "gap_lower_bound" => token(round_down(Q(b), big(10)^9)))
end

"""
Exact ground-energy enclosures a = (lo, hi) of E₀(n) and b of E₀(2n), and a cached
`logbar(L, β)` enclosing F_L(β) = log Tr exp[-β(H_L - a_L)] with a_L the lower end for ring L.
"""
function temperature_data(specs, n)
    a, b = Tuple(ground_energy_enclosure(specs[n])), Tuple(ground_energy_enclosure(specs[2n]))
    cache = Dict{Tuple{Int,Q},Interval{Float64}}()
    logbar(L, beta) = get!(cache, (L, beta)) do
        shifted_log_partition(specs[L], beta, L == n ? a[1] : b[1])
    end
    return (; a, b, logbar)
end

"""
Lower seed defects valid for every β in [lo, hi]. Spatial: F_L decreases (H_L - a_L ⪰ 0),
so log Pₛ ≤ F_2n(lo) - 2F_n(hi) - β(a_2n - 2a_n), evaluated over the whole interval; no
monotonicity of Pₛ is assumed. Thermal: d log P_τ/dβ = 2(⟨H⟩_β - ⟨H⟩_2β) ≥ 0, so P_τ ≤ P_τ(hi).
"""
function window_lower(data, n, lo, hi)
    beta = interval(Float64(lo), Float64(hi))
    @assert Q(Float64(lo)) == lo && Q(Float64(hi)) == hi
    logps = data.logbar(2n, lo) - 2data.logbar(n, hi) - beta * to_interval(data.b[1] - 2data.a[1])
    logpt = data.logbar(2n, 2hi) - 2data.logbar(2n, hi)
    return lower(-expm1(logps)), lower(-expm1(logpt))
end

"""
Reject every β ≥ B: with δ = lower E₀(2n) - 2·upper E₀(n) > 0,
log Pₛ(β) ≤ F_2n(β) - βδ ≤ F_2n(B) - Bδ.
"""
function tail_lower(data, n, B)
    delta = data.b[1] - 2data.a[2]
    delta > 0 || return nothing
    sl = lower(-expm1(data.logbar(2n, B) - to_interval(B * delta)))
    rejection(sl, big(0) // 1) === nothing && return nothing
    return Dict("beta_from" => token(B), "ground_mismatch_lower" => token(delta),
        "eps_s_lower" => token(sl), "eps_t_lower" => "0//1")
end

"""
    exclude_temperatures(specs, n) -> Dict

Prove that seed length `n` fails at every β > 0: find a tail start B, then bisect [0, B]
into windows each carrying a lower-defect rejection witness. Throws if unresolved.
"""
function exclude_temperatures(specs, n)
    data = temperature_data(specs, n)
    B = big(1) // 1
    tail = tail_lower(data, n, B)
    while tail === nothing && B < 1024
        B *= 2
        tail = tail_lower(data, n, B)
    end
    tail === nothing && error("could not exclude the infinite-temperature-parameter tail")
    windows = Any[]
    function split(lo, hi, depth)
        sl, tl = window_lower(data, n, lo, hi)
        witness = rejection(sl, tl)
        if witness !== nothing
            push!(windows, Dict("beta_from" => token(lo), "beta_to" => token(hi),
                "eps_s_lower" => token(sl), "eps_t_lower" => token(tl), "rejection" => witness))
        else
            depth < 20 || error("temperature window unresolved: $lo to $hi")
            split(lo, (lo + hi) / 2, depth + 1)
            split((lo + hi) / 2, hi, depth + 1)
        end
    end
    split(big(0) // 1, B, 0)
    return Dict("n0" => token(n), "windows" => windows, "tail" => tail)
end

"Check contiguous coverage of (0, ∞) and every lower-defect witness, including the tail."
function check_exclusion(specs, record)
    n = Int(q(record["n0"]))
    n in SEED_LENGTHS || error("invalid excluded seed length")
    data = temperature_data(specs, n)
    previous = big(0) // 1
    for row in record["windows"]
        lo, hi = q(row["beta_from"]), q(row["beta_to"])
        lo == previous && hi > lo || error("temperature coverage gap")
        sl, tl = window_lower(data, n, lo, hi)
        s, t = q(row["eps_s_lower"]), q(row["eps_t_lower"])
        s <= sl && t <= tl || error("temperature lower defect fails")
        rejection(s, t) == row["rejection"] || error("temperature rejection fails")
        previous = hi
    end
    tail = record["tail"]
    previous == q(tail["beta_from"]) || error("tail coverage gap")
    fresh = tail_lower(data, n, previous)
    fresh !== nothing || error("tail rejection fails")
    q(tail["ground_mismatch_lower"]) <= q(fresh["ground_mismatch_lower"]) ||
        error("ground mismatch fails")
    q(tail["eps_s_lower"]) <= q(fresh["eps_s_lower"]) || error("tail defect fails")
    q(tail["eps_t_lower"]) == 0 || error("invalid tail thermal defect")
    rejection(q(tail["eps_s_lower"]), big(0) // 1) !== nothing || error("tail witness fails")
    return true
end

"""
    certify_coupling(J) -> (record, certificate or nothing)

Classify every grid seed at coupling `J` from fresh verified spectra and, at the best
certified seed, produce the core gap certificate with `certify_gap` from that seed's
`seed_bounds` result (no spectrum is recomputed).
"""
function certify_coupling(J)
    start = time()
    check_frame(J)
    model = ladder_model(J)
    specs = verified_spectra(model)
    seeds = Dict((n, beta) => seed(model, specs, n, beta) for n in SEED_LENGTHS for beta in BETAS)
    rows = [classify(n, beta, seeds[(n, beta)]) for n in SEED_LENGTHS for beta in BETAS]
    best = best_row(rows)
    cert = nothing
    if best !== nothing
        n, beta = Int(q(best["n0"])), q(best["beta0"])
        cert = certify_gap(model; n0 = n, β0 = beta, rp = ladder_rp(model),
            seed = seeds[(n, beta)])
    end
    record = Dict{String,Any}("J" => token(J), "grid" => rows, "best" => best,
        "certificate" => cert === nothing ? nothing : basename(certificate_path(J)),
        "continuation" => cert === nothing ? nothing :
                          continuation(cert["seed"], CONTINUATION_DELTA),
        "all_beta_exclusions" => J in ALL_BETA_EXCLUDED ?
                                 [exclude_temperatures(specs, n) for n in SEED_LENGTHS] : nothing,
        "spectra" => Dict(string(L) => spectrum_summary(s) for (L, s) in specs),
        "seconds" => time() - start)
    return record, cert
end

"""
    produce() -> Dict

Certify every coupling in `COUPLINGS`: write one core certificate per certified coupling
(`results/certificate_J<tag>.json`, with provenance metadata) and the suite record
`results/suite.json`.
"""
function produce()
    start = time()
    records = Any[]
    meta = model_metadata(BASE, Dict("n0" => token.(SEED_LENGTHS),
        "beta_grid" => "1//4:1//4:16//1", "source_of_spectra" => "fresh complete verification");
        HASHED..., packages = (CertifiedGaps, IntervalArithmetic, JSON3), certified = true)
    for J in COUPLINGS
        record, cert = certify_coupling(J)
        if cert !== nothing
            cert["metadata"] = meta
            write_result(certificate_path(J), cert)
        end
        push!(records, record)
        println("COUPLING ", J, cert === nothing ? " no certified seed" :
            " gap >= $(cert["interval_bound_decimal"]) (interval), " *
            "> $(cert["gap_lower_bound_decimal"]) (exact chain)",
            " seconds=", round(record["seconds"]; digits = 1))
        flush(stdout)
    end
    suite = Dict("format" => FORMAT, "parameters" => records, "metadata" => meta,
        "validity" => "even L >= n0 rungs (N = 2L spins); one core certificate per coupling",
        "seconds" => time() - start)
    write_result(SUITE, suite; compact = true)
    return suite
end

"""
    check_decision(row, beta)

Re-derive one grid row's decision from its stored rational defects alone: a certified row
must reproduce its core interval gap bound, a rejected row must carry the exact `rejection`
witness of its lower defects. Needs no spectrum.
"""
function check_decision(row, beta)
    sl, tl = q(row["eps_s_lower"]), q(row["eps_t_lower"])
    su, tu = q(row["eps_s_upper"]), q(row["eps_t_upper"])
    0 <= sl <= su && 0 <= tl <= tu || error("defect bounds out of order")
    if row["status"] == "certified"
        b = certified_gap_bound(su, tu, beta).bound
        b !== nothing && q(row["interval_bound"]) == Q(b) || error("gap bound fails")
    elseif row["status"] == "rejected_by_recurrence"
        witness = rejection(sl, tl)
        witness !== nothing || error("rejection not proved")
        witness == row["rejection"] || error("rejection witness mismatch")
    else
        row["status"] == "unresolved" || error("unknown status")
    end
    return true
end

"Verify one stored grid row against fresh core defect enclosures, then its decision."
function check_row(row, b, beta)
    sl, tl = q(row["eps_s_lower"]), q(row["eps_t_lower"])
    su, tu = q(row["eps_s_upper"]), q(row["eps_t_upper"])
    sl <= max(0, Q(inf(b.defect_s))) && tl <= max(0, Q(inf(b.defect_t))) ||
        error("lower defects fail")
    Q(sup(b.defect_s)) <= su && Q(sup(b.defect_t)) <= tu || error("upper defects fail")
    return check_decision(row, beta)
end

"""
Verify one coupling's record and core certificate. `specs` are fresh verified spectra
(`check.jl --full`) or `nothing` (fast check: stored rational defects only).
"""
function check_coupling(J, record, specs)
    check_frame(J)
    model = ladder_model(J)
    verify_rp(model, ladder_rp(model)) || error("RP identity fails")
    expected = [(n, beta) for n in SEED_LENGTHS for beta in BETAS]
    length(record["grid"]) == length(expected) || error("incomplete grid")
    for (row, (n, beta)) in zip(record["grid"], expected)
        (Int(q(row["n0"])), q(row["beta0"])) == (n, beta) || error("grid mismatch")
        specs === nothing ? check_decision(row, beta) :
            check_row(row, seed(model, specs, n, beta), beta)
    end
    best = best_row(record["grid"])
    best == record["best"] || error("best seed mismatch")
    claim = CLAIMED_GAPS[J]
    if best === nothing
        record["certificate"] === nothing && record["continuation"] === nothing ||
            error("certificate without a certified seed")
        claim === nothing && all(r -> r["status"] == "rejected_by_recurrence", record["grid"]) ||
            error("claimed full grid rejection fails at J=$J")
    else
        record["certificate"] == basename(certificate_path(J)) ||
            error("certificate path mismatch")
        cert = parse_exact_json(read(certificate_path(J), String))
        verify_model_metadata(cert["metadata"], BASE; HASHED...)
        stored = decode_model(cert["model"])
        stored.bond == model.bond && stored.onsite == model.onsite &&
            stored.charge == model.charge && stored.flip == model.flip ||
            error("certificate model differs from ladder_model($J)")
        s = cert["seed"]
        (q(s["n0"]), q(s["beta0"]), q(s["eps_s"]), q(s["eps_t"])) ==
            (q(best["n0"]), q(best["beta0"]), q(best["eps_s_upper"]), q(best["eps_t_upper"])) ||
            error("certificate seed differs from the best grid seed")
        # The core checks the RP identity and the exact rational chain, and in the full mode
        # the seed defects from the fresh spectra of the (identical) decoded model.
        n = Int(q(s["n0"]))
        check_certificate(cert;
            spectra = specs === nothing ? false : Dict(n => specs[n], 2n => specs[2n]))
        continuation(cert["seed"], CONTINUATION_DELTA) == record["continuation"] ||
            error("continuation mismatch")
        # The README reports the core interval bound (certified_gap_bound of the seed); the
        # exact chain bound `gap_lower_bound` was re-derived by `check_certificate` above.
        b = certified_gap_bound(q(s["eps_s"]), q(s["eps_t"]), q(s["beta0"])).bound
        b !== nothing && Q(b) == Q(cert["interval_bound_decimal"]) ||
            error("interval bound mismatch at J=$J")
        claim !== nothing && Q(b) >= claim || error("claimed gap bound fails at J=$J")
    end
    exclusions = record["all_beta_exclusions"]
    if J in ALL_BETA_EXCLUDED
        q.(getindex.(exclusions, "n0")) == collect(SEED_LENGTHS) ||
            error("exclusion lengths mismatch")
        specs === nothing || foreach(r -> check_exclusion(specs, r), exclusions)
    else
        exclusions === nothing || error("unexpected exclusion")
    end
    return true
end

"""
    check(; full = false) -> true

Verify `results/suite.json` and every core certificate it names. The fast mode (`check.jl`)
checks hashes, the exact frame and RP identities, every grid decision from its stored rational
defects, the best seeds, each core certificate with `check_certificate(cert; spectra = false)`
(RP and the exact rational chain), the continuations and the README claims. `full = true`
(`check.jl --full`) also recomputes every verified spectrum, every grid defect enclosure, each
certificate's seed through `check_certificate` with those fresh spectra, and the all-β
exclusions. Throws on any mismatch.
"""
function check(; full = false)
    start = time()
    suite = parse_exact_json(read(SUITE, String))
    suite["format"] == FORMAT || error("unknown suite format")
    verify_model_metadata(suite["metadata"], BASE; HASHED...)
    q.(getindex.(suite["parameters"], "J")) == Q.(COUPLINGS) || error("couplings mismatch")
    certified = [J for J in COUPLINGS if CLAIMED_GAPS[J] !== nothing]
    on_disk = filter(f -> startswith(f, "certificate"), readdir(joinpath(BASE, "results")))
    Set(on_disk) == Set(basename.(certificate_path.(certified))) ||
        error("unexpected or missing certificate files")
    for record in suite["parameters"]
        J = Rational{Int}(q(record["J"]))
        specs = full ? verified_spectra(ladder_model(J)) : nothing
        check_coupling(J, record, specs)
        full && println("checked J = $J; ", round(time() - start; digits = 1), " s")
    end
    println(full ? "PASS (full, fresh spectra): " : "PASS (fast, no spectra): ",
        "hashes, frame, RP, grid decisions, core certificates, exclusions, continuations ",
        "and README claims; ", round(time() - start; digits = 2), " seconds")
    return true
end
