# The certificate suite: fresh verified spectra at every coupling, every grid seed classified
# rigorously, all-beta exclusions where claimed, local coupling continuation, and a checker
# that recomputes all of it. No numerical scan or older certificate is an input.

"Exact multiple used for outward-rounded seed defects."
const GRID = big(10)^8
"Seed temperatures β₀ = k/4, k = 1, …, 64."
const BETAS = [big(k) // 4 for k in 1:64]
"Even seed lengths reachable with complete spectra of rings of at most 8 rungs."
const SEED_LENGTHS = (2, 4)
"""
Couplings in the suite: the four certified ones, the two all-β exclusions, and the two
endpoints of the finite-grid threshold bracket (proposed by `scripts/threshold_bisection.jl`).
"""
const COUPLINGS = [4 // 1, 7 // 2, 3 // 1, 5 // 2, 2 // 1, 3 // 2,
    326133 // 131072, 163069 // 65536]
"Couplings at which seeds n₀ ∈ {2, 4} are excluded for every β > 0."
const ALL_BETA_EXCLUDED = (2 // 1, 3 // 2)
"""
README headline claims checked by `check`: a certified gap lower bound, or `nothing` where
every grid seed must be rigorously rejected.
"""
const CLAIMED_GAPS = Dict{Rational{Int},Union{Q,Nothing}}(
    4 // 1 => 1102403091 // 1000000000, 7 // 2 => 77864819 // 100000000,
    3 // 1 => 227094667 // 500000000, 5 // 2 => 42433721 // 1000000000,
    2 // 1 => nothing, 3 // 2 => nothing, 326133 // 131072 => nothing,
    163069 // 65536 => 306233 // 200000000)
"Coupling half-width of the continuation intervals around certified couplings."
const CONTINUATION_DELTA = big(1) // 10000
const FORMAT = "heisenberg-ladder verified suite v1"
const VALIDITY = "even L >= n0 rungs; N=2L spin sites, N divisible by 4"
"The certificate file produced by `certify.jl` and verified by `check.jl`."
const CERTIFICATE = joinpath(BASE, "results", "certificate.json")

q(x::AbstractString) = parse(Q, x)
q(x::Integer) = Q(x)
token(x) = string(Q(x))
endpoints(x) = [token(inf(x)), token(sup(x))]

"""
    provenance(config; certified = false)

Run metadata with SHA-256 of the shared core sources, this model's library, its entry points
and environment. `check` re-hashes exactly these files. Exploratory scripts and tests are
deliberately not hashed, so editing them does not invalidate the certificate.
"""
function provenance(config; certified = false)
    model = "models/heisenberg_ladder"
    hashes = source_hashes(ROOT; dirs = ("julia/src", "$model/src"),
        files = ("julia/Project.toml", "$model/Project.toml", "$model/Manifest.toml",
            "$model/certify.jl", "$model/check.jl", "notes/thermal-mpo-gap-scope.tex"))
    return run_metadata(config; source_sha256 = hashes,
        packages = (CertifiedGaps, IntervalArithmetic, JSON3),
        status = certified ? "certified" : "numerical, not certified", certified)
end

"Write `data` as single-line JSON, atomically (the pretty form is ~2.3 times larger)."
function write_compact(path, data)
    mkpath(dirname(abspath(path)))
    temporary = path * ".tmp"
    open(io -> (JSON3.write(io, data); println(io)), temporary, "w")
    mv(temporary, path; force = true)
    return path
end

"Enclose both true seed defects using every verified eigenvalue and its multiplicity."
function defects(specs, n, beta)
    a = enclosed_logZ(specs[n], beta)
    b = enclosed_logZ(specs[2n], beta)
    c = enclosed_logZ(specs[2n], 2beta)
    return -expm1(b - 2a), -expm1(c - 2b), [a, b, c]
end

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

"""
    classify(specs, n, beta) -> Dict

Classify one seed: `certified` (upper defects close the core recurrence),
`rejected_by_recurrence` (lower defects carry a `rejection` witness) or `unresolved`.
A failed upper-bound attempt is never called a rejection.
"""
function classify(specs, n, beta)
    s, t, logs = defects(specs, n, beta)
    su = max(1 // GRID, round_up(Q(sup(s)), GRID))
    tu = max(1 // GRID, round_up(Q(sup(t)), GRID))
    sl = max(big(0) // 1, round_down(Q(inf(s)), GRID))
    tl = max(big(0) // 1, round_down(Q(inf(t)), GRID))
    answer = su < 1 && tu < 1 ? certified_gap_bound(su, tu, beta).bound : nothing
    witness = rejection(sl, tl)
    answer !== nothing && witness !== nothing && error("inconsistent classification")
    status = answer !== nothing ? "certified" :
        witness !== nothing ? "rejected_by_recurrence" : "unresolved"
    record = Dict{String,Any}("n0" => token(n), "beta0" => token(beta), "status" => status,
        "eps_s_upper" => token(su), "eps_t_upper" => token(tu),
        "eps_s_lower" => token(sl), "eps_t_lower" => token(tl))
    if answer !== nothing
        record["gap_lower_bound"] = token(round_down(Q(answer), big(10)^9))
        record["interval_bound_exact"] = token(answer)
        record["logZ"] = endpoints.(logs)
    end
    witness !== nothing && (record["rejection"] = witness)
    return record
end

"""
    continuation(seed, delta) -> Dict or nothing

Enlarge a certified seed to every coupling within `delta`. The centred rung sum satisfies
‖W_L + L/4‖ ≤ L/2, so each log Z_L moves by at most βδL/2 after a scalar shift that cancels in
both purities; this gives s_δ = 1-(1-s)e^{-2βnδ} and t_δ = 1-(1-t)e^{-4βnδ}.
"""
function continuation(seed, delta)
    n, beta = Int(q(seed["n0"])), q(seed["beta0"])
    s, t = q(seed["eps_s_upper"]), q(seed["eps_t_upper"])
    su = round_up(Q(sup(1 - (1 - to_interval(s)) * exp(-to_interval(2beta * n * delta)))), GRID)
    tu = round_up(Q(sup(1 - (1 - to_interval(t)) * exp(-to_interval(4beta * n * delta)))), GRID)
    b = certified_gap_bound(su, tu, beta).bound
    b === nothing && return nothing
    return Dict("delta" => token(delta), "eps_s_upper" => token(su), "eps_t_upper" => token(tu),
        "gap_lower_bound" => token(round_down(Q(b), big(10)^9)))
end

"Ground energy bounds from the complete spectrum, with no ordering assumption on radii."
ground_bounds(s) = (minimum(Q(e) - Q(r) for (e, r) in zip(s.energies, s.radii)),
    minimum(Q(e) + Q(r) for (e, r) in zip(s.energies, s.radii)))

"""
Ground-energy bounds for rings n and 2n, and a cached `logbar(L, β)` enclosing
F_L(β) = log Tr exp[-β(H_L - a_L)] with a_L a lower bound on E₀(L).
"""
function temperature_data(specs, n)
    a, b = ground_bounds(specs[n]), ground_bounds(specs[2n])
    cache = Dict{Tuple{Int,Q},Interval{Float64}}()
    function logbar(L, beta)
        get!(cache, (L, beta)) do
            shift = L == n ? a[1] : b[1]
            enclosed_logZ(specs[L], beta) + to_interval(beta * shift)
        end
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
    sl = max(big(0) // 1, round_down(Q(inf(-expm1(logps))), GRID))
    tl = max(big(0) // 1, round_down(Q(inf(-expm1(logpt))), GRID))
    return sl, tl
end

"""
Reject every β ≥ B: with δ = lower E₀(2n) - 2·upper E₀(n) > 0,
log Pₛ(β) ≤ F_2n(β) - βδ ≤ F_2n(B) - Bδ.
"""
function tail_lower(data, n, B)
    delta = data.b[1] - 2data.a[2]
    delta > 0 || return nothing
    logps = data.logbar(2n, B) - to_interval(B * delta)
    sl = max(big(0) // 1, round_down(Q(inf(-expm1(logps))), GRID))
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

"Fresh verified spectra of the rings of 2, 4 and 8 rungs."
verified_spectra(J) = Dict(L => full_spectrum(L, J; verified = true) for L in (2, 4, 8))

"Highest certified grid row by gap lower bound, or `nothing`."
function best_row(rows)
    good = filter(r -> r["status"] == "certified", rows)
    return isempty(good) ? nothing : argmax(r -> q(r["gap_lower_bound"]), good)
end

"Recompute one coupling; the output includes every beta-grid decision, including failures."
function certify_coupling(J)
    start = time()
    @assert exact_rp(J)
    model, rp = rp_model(J)
    specs = verified_spectra(J)
    rows = [classify(specs, n, beta) for n in SEED_LENGTHS for beta in BETAS]
    best = best_row(rows)
    cont = best === nothing ? nothing : continuation(best, CONTINUATION_DELTA)
    excluded = J in ALL_BETA_EXCLUDED ?
        [exclude_temperatures(specs, n) for n in SEED_LENGTHS] : nothing
    summaries = Dict(string(L) => Dict("dimension" => token(s.dimension),
        "max_radius" => token(round_up(Q(maximum(s.radii)), GRID)), "seconds" => s.seconds,
        "block_count" => token(length(s.rows))) for (L, s) in specs)
    return Dict("J" => token(J), "model" => encode_model(model), "rp" => encode_rp(rp),
        "grid" => rows, "best" => best, "continuation" => cont, "spectra" => summaries,
        "all_beta_exclusions" => excluded, "seconds" => time() - start)
end

"""
    produce(path = CERTIFICATE) -> Dict

Certify every coupling in `COUPLINGS` from fresh verified spectra and write the certificate.
"""
function produce(path = CERTIFICATE)
    start = time()
    records = Any[]
    for J in COUPLINGS
        result = certify_coupling(J)
        push!(records, result)
        println("COUPLING ", J, " best=", result["best"] === nothing ? "none" :
            result["best"]["gap_lower_bound"], " seconds=", result["seconds"])
        flush(stdout)
    end
    cert = Dict("format" => FORMAT, "parameters" => records,
        "metadata" => provenance(Dict("n0" => token.(SEED_LENGTHS),
            "beta_grid" => "1//4:1//4:16//1",
            "source_of_spectra" => "fresh complete verification"); certified = true),
        "validity" => VALIDITY, "seconds" => time() - start)
    write_compact(path, cert)
    return cert
end

"""
    check_decision(row, beta)

Re-derive one grid row's decision from its stored rational defects alone: a certified row must
close the core recurrence with its upper defects and its stored gap, a rejected row must carry
the exact `rejection` witness of its lower defects. Needs no spectrum.
"""
function check_decision(row, beta)
    sl, tl = q(row["eps_s_lower"]), q(row["eps_t_lower"])
    su, tu = q(row["eps_s_upper"]), q(row["eps_t_upper"])
    0 <= sl <= su && 0 <= tl <= tu || error("defect bounds out of order")
    if row["status"] == "certified"
        b = certified_gap_bound(su, tu, beta).bound
        b !== nothing && 0 < q(row["gap_lower_bound"]) <= Q(b) || error("gap bound fails")
        q(row["interval_bound_exact"]) <= Q(b) || error("interval bound mismatch")
    elseif row["status"] == "rejected_by_recurrence"
        witness = rejection(sl, tl)
        witness !== nothing || error("rejection not proved")
        witness == row["rejection"] || error("rejection witness mismatch")
    else
        row["status"] == "unresolved" || error("unknown status")
    end
    return true
end

"Verify one stored grid row against fresh defect enclosures."
function check_row(row, specs, n, beta)
    (Int(q(row["n0"])), q(row["beta0"])) == (n, beta) || error("grid mismatch")
    s, t, logs = defects(specs, n, beta)
    sl, tl = q(row["eps_s_lower"]), q(row["eps_t_lower"])
    su, tu = q(row["eps_s_upper"]), q(row["eps_t_upper"])
    sl <= max(0, Q(inf(s))) && tl <= max(0, Q(inf(t))) || error("lower defects fail")
    Q(sup(s)) <= su && Q(sup(t)) <= tu || error("upper defects fail")
    check_decision(row, beta)
    if row["status"] == "certified"
        for (saved, fresh) in zip(row["logZ"], logs)
            # Independent eigensolves need only overlap; the seeds are checked above.
            q(saved[1]) <= Q(sup(fresh)) && Q(inf(fresh)) <= q(saved[2]) ||
                error("log Z enclosure mismatch")
        end
    end
end

"Verify the README headline for one coupling against the (already checked) grid."
function check_claim(J, result)
    claim = CLAIMED_GAPS[J]
    if claim === nothing
        all(r -> r["status"] == "rejected_by_recurrence", result["grid"]) ||
            error("claimed full grid rejection fails at J=$J")
    else
        result["best"] !== nothing && q(result["best"]["gap_lower_bound"]) >= claim ||
            error("claimed gap bound fails at J=$J")
    end
    return true
end

"Format, validity claim, source hashes and coupling list of a parsed certificate."
function check_header(cert)
    cert["format"] == FORMAT || error("unknown certificate format")
    cert["validity"] == VALIDITY || error("incorrect validity claim")
    for (file, hash) in cert["metadata"]["source_sha256"]
        full = joinpath(ROOT, file)
        isfile(full) || error("hashed source missing: $file")
        sha256_file(full) == hash || error("source mismatch: $file")
    end
    q.(getindex.(cert["parameters"], "J")) == Q.(COUPLINGS) || error("couplings mismatch")
    return true
end

"Exact RP of coupling `J` and agreement of the stored model and RP factors with it."
function check_model(J, result)
    @assert exact_rp(J)
    model, rp = rp_model(J)
    stored = decode_model(result["model"])
    stored.d == model.d && stored.bond == model.bond && stored.onsite == model.onsite ||
        error("model mismatch")
    verify_rp(stored, decode_rp(result["rp"])) || error("RP identity fails")
    return true
end

"""
    check_stored(path = CERTIFICATE) -> true

Fast check (`check.jl`, seconds): everything `check` verifies that needs no spectrum. It
re-derives the source hashes, the exact RP identities and stored model, every grid decision
from its stored rational defects, the best seed, the coupling continuation and the README
claims. The defect enclosures themselves, the log Z intervals and the all-β exclusions are
re-established only by `check` (`check.jl --full`). Throws on any mismatch.
"""
function check_stored(path = CERTIFICATE)
    start = time()
    cert = parse_exact_json(read(path, String))
    check_header(cert)
    for result in cert["parameters"]
        J = Rational{Int}(q(result["J"]))
        check_model(J, result)
        expected = [(n, beta) for n in SEED_LENGTHS for beta in BETAS]
        length(result["grid"]) == length(expected) || error("incomplete grid")
        for (row, (n, beta)) in zip(result["grid"], expected)
            (Int(q(row["n0"])), q(row["beta0"])) == (n, beta) || error("grid mismatch")
            check_decision(row, beta)
        end
        best = best_row(result["grid"])
        best == result["best"] || error("best seed mismatch")
        if result["continuation"] !== nothing
            fresh = continuation(best, q(result["continuation"]["delta"]))
            fresh == result["continuation"] || error("continuation mismatch")
        end
        exclusions = result["all_beta_exclusions"]
        if J in ALL_BETA_EXCLUDED
            q.(getindex.(exclusions, "n0")) == collect(SEED_LENGTHS) ||
                error("exclusion lengths mismatch")
        else
            exclusions === nothing || error("unexpected exclusion")
        end
        check_claim(J, result)
    end
    println("PASS (fast, no spectra): hashes, RP, grid decisions, best seeds, continuation ",
        "and README claims; ", round(time() - start; digits = 2), " seconds")
    return true
end

"""
    check(path = CERTIFICATE) -> true

Verify the certificate from scratch (`check.jl --full`): source hashes, exact RP, fresh
complete spectra, every seed inequality and grid decision, the core recurrence, all-β
exclusions, coupling continuation and the README claims. Throws on any mismatch.
"""
function check(path = CERTIFICATE)
    start = time()
    cert = parse_exact_json(read(path, String))
    check_header(cert)
    for result in cert["parameters"]
        J = Rational{Int}(q(result["J"]))
        check_model(J, result)
        specs = verified_spectra(J)
        expected = [(n, beta) for n in SEED_LENGTHS for beta in BETAS]
        length(result["grid"]) == length(expected) || error("incomplete grid")
        for (row, (n, beta)) in zip(result["grid"], expected)
            check_row(row, specs, n, beta)
        end
        best = best_row(result["grid"])
        best == result["best"] || error("best seed mismatch")
        if result["continuation"] !== nothing
            fresh = continuation(best, q(result["continuation"]["delta"]))
            fresh == result["continuation"] || error("continuation mismatch")
        end
        if J in ALL_BETA_EXCLUDED
            exclusions = result["all_beta_exclusions"]
            q.(getindex.(exclusions, "n0")) == collect(SEED_LENGTHS) ||
                error("exclusion lengths mismatch")
            foreach(r -> check_exclusion(specs, r), exclusions)
        else
            result["all_beta_exclusions"] === nothing || error("unexpected exclusion")
        end
        for (L, s) in specs
            summary = result["spectra"][string(L)]
            Int(q(summary["dimension"])) == s.dimension == 4^L || error("dimension mismatch")
            Int(q(summary["block_count"])) == length(s.rows) || error("block count mismatch")
            Q(maximum(s.radii)) <= q(summary["max_radius"]) || error("spectral radius mismatch")
        end
        check_claim(J, result)
    end
    println("PASS: fresh RP, complete spectra, seed bounds, grid decisions, exclusions, ",
        "continuation and README claims; ", round(time() - start; digits = 2), " seconds")
    return true
end
