# The certificate suite on top of the core pipeline. Per coupling: verified spectra of the
# rings of 2, 4 and 8 rungs (`verified_ring_spectrum`), every grid seed classified from the
# core seed defects (`seed_bounds`), a core gap certificate (`certify_point`) at the best
# certified seed, the core all-β exclusions (`Exclusions.exclude_all_beta`) where no seed
# exists, and the model-specific coupling continuation and README claims. No numerical scan or
# older result is an input.

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
    4 // 1 => big(1102404333) // 10^9, 7 // 2 => big(778648931) // 10^9,
    3 // 1 => big(454189767) // 10^9, 5 // 2 => big(42434526) // 10^9,
    2 // 1 => nothing, 3 // 2 => nothing, 326133 // 131072 => nothing,
    163069 // 65536 => big(1540056) // 10^9)
"Coupling half-width of the continuation intervals around certified couplings."
const CONTINUATION_DELTA = big(1) // 10000
"Grid of the lower defect bounds and continuation defects (the core seeds also use 10⁻¹⁵)."
const GRID = big(10)^15
const FORMAT = "heisenberg-ladder suite v3 (core certificates and exclusions)"
"Model-specific suite record: grid decisions, exclusions, continuations and summaries."
const SUITE = joinpath(BASE, "results", "suite.json")
"The Lean gap theorem of the ladder, whose `points` hardcode seed defects."
const LEAN_FILE = joinpath(BASE, "..", "..", "lean", "CertifiedGaps", "Models",
    "HeisenbergLadder.lean")

q(x::AbstractString) = parse(Q, x)
q(x::Integer) = Q(x)
token(x) = string(Q(x))
"File tag of a coupling, e.g. 7//2 -> \"7_2\"."
coupling_tag(J) = denominator(J) == 1 ? string(numerator(J)) :
                  string(numerator(J), "_", denominator(J))
"Core certificate of a certified coupling."
certificate_path(J) = joinpath(BASE, "results", "certificate_J$(coupling_tag(J)).json")
"Exact parameters of the coupling `J`, recorded in its certificate."
point_params(J) = Dict("J" => J)

"""
    lean_seed(J, beta) -> (; s, t) or nothing

The seed defects that the Lean theorem (`LEAN_FILE`, read by `lean_points`) assumes at the
coupling `J`, or `nothing` if it has no such point. Its β₀ must be the certificate's `beta`.
"""
function lean_seed(J, beta)
    for p in lean_points(LEAN_FILE)
        p.params[1:end-1] == [J] || continue
        p.params[end] == beta || error("the Lean point J = $J has β₀ = $(p.params[end])")
        return (; p.s, p.t)
    end
    return nothing
end

"Fresh verified spectra (core `BlockSpectrum`s) of the rings of 2, 4 and 8 rungs."
verified_spectra(model) = Dict(L => verified_ring_spectrum(model, L) for L in RINGS)

"Core seed defects of the seed (n, β), reusing the precomputed spectra of rings n and 2n."
seed(model, specs, n, beta) =
    seed_bounds(model, n, beta; spectra = Dict(n => specs[n], 2n => specs[2n]))

"Lower end of a defect enclosure, rounded down to `GRID` and clipped at zero."
lower(x::Interval) = max(big(0) // 1, round_down(Q(inf(x)), GRID))

"""
    classify(n, beta, seed) -> Dict

Classify one grid seed from its core `seed_bounds` result: `certified` (the upper defects
close the core recurrence), `rejected_by_recurrence` (the lower defects carry an
`Exclusions.rejection` witness) or `unresolved`. A failed upper-bound attempt is never called
a rejection.
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
The centered rung sum satisfies ‖W_L + L/4‖ ≤ L/2, so each log Z_L moves by at most βδL/2 after
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
    certify_coupling(J) -> (record, certificate or nothing)

Classify every grid seed at coupling `J` from fresh verified spectra and, at the best
certified seed, write the core certificate with `certify_point` from that seed's
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
        cert = certify_point(certificate_path(J), model; n0 = n, β0 = beta,
            params = point_params(J), model_dir = BASE, rp = ladder_rp(model),
            seed = seeds[(n, beta)], packages = (IntervalArithmetic, JSON3))
    end
    record = Dict{String,Any}("J" => token(J), "grid" => rows, "best" => best,
        "certificate" => cert === nothing ? nothing : basename(certificate_path(J)),
        "continuation" => cert === nothing ? nothing :
                          continuation(cert["seed"], CONTINUATION_DELTA),
        "all_beta_exclusions" => J in ALL_BETA_EXCLUDED ?
            [exclude_all_beta(specs[n], specs[2n], n) for n in SEED_LENGTHS] : nothing,
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
        packages = (CertifiedGaps, IntervalArithmetic, JSON3), certified = true)
    for J in COUPLINGS
        record, cert = certify_coupling(J)
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
        path = certificate_path(J)
        record["certificate"] == basename(path) || error("certificate path mismatch")
        # The core checks the hashes, the point's model, RP data, seed, params and Lean seed,
        # the RP identity, the exact rational chain and the interval bound, and in the full
        # mode the seed defects from the fresh spectra of the identical model.
        n, beta = Int(q(best["n0"])), q(best["beta0"])
        check_point(path, model; n0 = n, β0 = beta, params = point_params(J), model_dir = BASE,
            rp = ladder_rp(model), lean = lean_seed(J, beta),
            spectra = specs === nothing ? false : Dict(n => specs[n], 2n => specs[2n]))
        cert = read_result(path)
        s = cert["seed"]
        (q(s["eps_s"]), q(s["eps_t"])) == (q(best["eps_s_upper"]), q(best["eps_t_upper"])) ||
            error("certificate seed differs from the best grid seed")
        continuation(s, CONTINUATION_DELTA) == record["continuation"] ||
            error("continuation mismatch")
        claim !== nothing && Q(cert["interval_bound_decimal"]) >= claim ||
            error("claimed gap bound fails at J=$J")
    end
    exclusions = record["all_beta_exclusions"]
    if J in ALL_BETA_EXCLUDED
        q.(getindex.(exclusions, "n0")) == collect(SEED_LENGTHS) ||
            error("exclusion lengths mismatch")
        for r in exclusions
            n = Int(q(r["n0"]))
            specs === nothing ? check_exclusion(r) :
                check_exclusion(r; small = specs[n], large = specs[2n])
        end
    else
        exclusions === nothing || error("unexpected exclusion")
    end
    return true
end

"""
    check(; full = false) -> true

Verify `results/suite.json` and every core certificate it names. The fast mode (`check.jl`)
checks hashes, the exact frame and RP identities, every grid decision from its stored rational
defects, the best seeds, each core certificate with `check_point(...; spectra = false)` (the
point's definition, the Lean seed, RP, the exact rational chain and the interval bound), the
continuations, the coverage and witnesses of the all-β exclusions and the README claims.
`full = true` (`check.jl --full`) also recomputes every verified spectrum, every grid defect
enclosure, each certificate's seed through `check_point` with those fresh spectra, and every
exclusion bound. Throws on any mismatch.
"""
function check(; full = false)
    start = time()
    suite = parse_exact_json(read(SUITE, String))
    suite["format"] == FORMAT || error("unknown suite format")
    verify_model_metadata(suite["metadata"], BASE)
    q.(getindex.(suite["parameters"], "J")) == Q.(COUPLINGS) || error("couplings mismatch")
    certified = [J for J in COUPLINGS if CLAIMED_GAPS[J] !== nothing]
    all(p -> length(p.params) == 2 && p.params[1] in certified, lean_points(LEAN_FILE)) ||
        error("a Lean point is not a certified coupling")
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
        "hashes, frame, RP, grid decisions, core certificates, Lean seed, exclusions, ",
        "continuations and README claims; ", round(time() - start; digits = 2), " seconds")
    return true
end
