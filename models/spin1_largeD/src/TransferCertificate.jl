"""
    TransferCertificate

Stage 2 of the spin-1 large-D certificate: rigorous seeds beyond exact diagonalisation.

From the complete, verified L = 6 and L = 12 spectral enclosures stored in `results/`, this
module encloses log Z, lifts the seed length by transfer-operator squaring (n₀ = 6·2ʲ),
transports seeds to exact D intervals by a centred Weyl bound, and passes the exact rational
defects to `CertifiedGaps.certified_gap_bound`. It also builds and replays the finite L = 8, 10
bridge. The inequalities are proved in `METHOD.md`.
"""
module TransferCertificate

# model.jl stays at the top level: the stage-1 certificates record its SHA-256 under that path.
include("../model.jl")
include("spectra.jl")

const ROOT = dirname(@__DIR__)
const Q = Rational{BigInt}
const box = CertifiedGaps.to_interval
const FORMAT = "spin1_largeD transfer and Weyl certificate v1"
const VALIDITY = "unique ground state and gap >= bound for every even L >= n0 sites"
const KEEP = 8192
"The case ids `scripts/transfer_certify.jl` writes; `check_certificate` requires exactly these."
const EXPECTED_CASES = ["reference-D4//1", "survey-D4-beta3",
    ["space-lift-$k-D4//1" for k in 1:4]..., "transport-threshold-D4//1",
    "reference-D1603//500", "transport-threshold-D1603//500"]
exact(x) = string(Q(x))
endpoints(x) = [exact(inf(x)), exact(sup(x))]
upper(x) = box(max(0, sup(x)))

"File stem of the stored results for an exact D, e.g. `D1603_500`."
stem(D) = "D$(numerator(D))_$(denominator(D))"

"""
    reference_spectrum(D, L)

Load the complete verified spectrum (centres `E`, radii `R`) at ring length `L`, after
checking its bytes against the SHA-256 recorded in the stage-1 certificate.
"""
function reference_spectrum(D, L)
    cert = parse_exact_json(read(joinpath(ROOT, "results", "certificate_$(stem(D)).json"),
        String))
    name = "enclosed_$(stem(D))_L$L"
    text = enclosure_text(joinpath(ROOT, "results"), name)
    hash = text_sha256(text)
    hash == cert["enclosure_hashes"]["results/$name.json"] || error("spectrum hash mismatch")
    data = JSON3.read(text, Dict{String,Any})
    data["rigorous"] && data["complete"] || error("incomplete/unverified spectrum")
    parse(Q, data["D"]) == D && data["L"] == L || error("wrong model parameters")
    E, R = Float64.(data["energies"]), Float64.(data["radii"])
    length(E) == length(R) == 3^L || error("wrong spectral dimension")
    sum(s["dimension"] * s["multiplicity"] for s in data["sectors"]) == 3^L ||
        error("wrong sector multiplicities")
    all(isfinite, E) && all(r -> isfinite(r) && r >= 0, R) && issorted(E) ||
        error("invalid spectral entries")
    return (; D, L, E, R, hash)
end

"""
    logz(spec, beta; keep = $KEEP)

Interval enclosure of log Z(β) from a complete spectrum: the lowest `keep` levels are summed
with interval exponentials, the rest by an explicit tail [0, (N−K) e^{−β(E_{K+1}−r_max)}].
"""
function logz(spec, beta; keep = KEEP)
    beta > 0 || throw(ArgumentError("positive beta required"))
    E, R = spec.E, spec.R
    k = min(keep, length(E))
    k > 0 || throw(ArgumentError("at least one level required"))
    b, shift = box(beta), box(E[1])
    total = interval(0.0)
    for j in 1:k
        total += exp(-b * (box(E[j]) + interval(-R[j], R[j]) - shift))
    end
    if k < length(E)
        # Centres are sorted, so every omitted true level is >= E[k+1] - max omitted radius.
        floor = box(E[k+1]) - box(maximum(@view R[k+1:end]))
        tail = (length(E) - k) * exp(-b * (floor - shift))
        total += interval(0.0, sup(tail))
    end
    return log(total) - b * shift
end

"""
    transported_logz(z, L, beta, Dref, lo, hi)

Centred Weyl transport: enclosure of log Z for H(D) − L(D−Dref)/2 at every D in [lo, hi],
from an enclosure `z` at Dref. The constant shift cancels in both purities.
"""
function transported_logz(z, L, beta, Dref, lo, hi)
    lo <= hi || throw(ArgumentError("reversed D interval"))
    radius = sup(box(L) * box(beta) * box(max(abs(lo-Dref), abs(hi-Dref))) / 2)
    return z + interval(-radius, radius)
end

"The purity-squaring bound s ↦ s²/(2(1−s)²), capped by the trivial defect bound one."
function square_defect(s)
    0 <= inf(s) <= sup(s) <= 1 || throw(ArgumentError("invalid defect"))
    sup(s) == 1 && return interval(1.0)
    return box(min(1, sup(s^2 / (2 * (1-s)^2))))
end

"""
    seed_ladder(z6, z12, z6twice, z12twice; max_lifts = 5)

Seeds at n₀ = 6·2ʲ, j = 0…max_lifts, from log Z₆(β), log Z₁₂(β), log Z₆(2β), log Z₁₂(2β).
Each lift squares both spatial defects and updates the thermal one by the exact coupling
identity P_τ(4n,β) = P_s(2n,2β) P_τ(2n,β)² / P_s(2n,β)² with P_s ≤ 1.
"""
function seed_ladder(z6, z12, z6twice, z12twice; max_lifts = 5)
    s = upper(1 - exp(z12 - 2z6))
    s2 = upper(1 - exp(z12twice - 2z6twice))
    t = upper(1 - exp(z12twice - 2z12))
    rows = [(; n0 = 6, eps_s = Q(sup(s)), eps_t = Q(sup(t)), lifts = 0)]
    for j in 1:max_lifts
        s, s2 = square_defect(s), square_defect(s2)
        t = upper(1 - (1-s2) * (1-t)^2)
        push!(rows, (; n0 = 6 * 2^j, eps_s = Q(sup(s)), eps_t = Q(sup(t)), lifts = j))
    end
    return rows
end

"""
    evaluate(a, b, beta; lo = a.D, hi = lo, max_lifts = 5, keep = $KEEP)

All seeds of the ladder for reference spectra `a` (L = 6) and `b` (L = 12), transported to
the exact interval [lo, hi]. See `evaluate_logs`.
"""
function evaluate(a, b, beta; lo = a.D, hi = lo, max_lifts = 5, keep = KEEP)
    a.L == 6 && b.L == 12 && a.D == b.D || error("expected matching L=6,12 spectra")
    rp_identity(lo) && rp_identity(hi) || error("RP identity failed")
    base = [logz(a, beta; keep), logz(b, beta; keep),
        logz(a, 2beta; keep), logz(b, 2beta; keep)]
    return evaluate_logs(base, beta, a.D, lo, hi; max_lifts)
end

"""
    evaluate_logs(base, beta, Dref, lo, hi; max_lifts = 5)

Seed rows from four cached log Z enclosures; `result.bound === nothing` marks a seed outside
the admissible region (a failure of these estimates, not a statement about the gap).
"""
function evaluate_logs(base, beta, Dref, lo, hi; max_lifts = 5)
    zs = [transported_logz(z, L, temp, Dref, lo, hi) for (z, L, temp) in
        zip(base, (6, 12, 6, 12), (beta, beta, 2beta, 2beta))]
    rows = seed_ladder(zs...; max_lifts)
    failed = (; bound = nothing, upper = nothing, j_admissible = nothing, log_rJ = nothing)
    # Exact pre-screen: the core recurrence needs combined thermal purity above one half.
    admissible(row) = row.eps_s < 1 && row.eps_t < 1 &&
        (1-row.eps_s)^2 * (1-row.eps_t) > 1//2
    return [(; row..., beta, lo, hi, logZ = zs, result = admissible(row) ?
        certified_gap_bound(row.eps_s, row.eps_t, beta) : failed) for row in rows]
end

"""
    finite_gap(spec, lo, hi)

Enclosure of the ring's spectral gap at every D in [lo, hi], including the Weyl shift L|δD|.
"""
function finite_gap(spec, lo, hi)
    width = box(spec.L) * box(max(abs(lo-spec.D), abs(hi-spec.D)))
    # Blockwise enclosures can cross after sorting the combined list. A common radius
    # encloses the globally sorted eigenvalues without assuming their block identities.
    r = 2 * box(maximum(spec.R)) + width
    return box(spec.E[2]) - box(spec.E[1]) + interval(-sup(r), sup(r))
end

"""
    proof_hashes()

SHA-256 of the model-local files that define the stage-2 proof and its environment (keyed
relative to this model directory). The checker requires exactly these keys and digests.
"""
function proof_hashes()
    return source_hashes(ROOT; dirs = ("src",),
        files = ("model.jl", "scripts/transfer_certify.jl", "check.jl", "Project.toml",
            "Manifest.toml"))
end

"""
    core_hashes()

SHA-256 of the shared core (`julia/src`, `julia/Project.toml`, `julia/Manifest.toml`, keyed by
repository path), stored as `run_metadata["core_sha256"]`. Provenance only: the checker does
not compare it with the checkout, since it reruns the core arithmetic itself.
"""
function core_hashes()
    repo = normpath(joinpath(ROOT, "..", ".."))
    return source_hashes(repo; dirs = ("julia/src",),
        files = ("julia/Project.toml", "julia/Manifest.toml"))
end

"Byte hashes of the stage-1 spectral evidence and the certificates that anchor it."
function input_hashes()
    names = String[]
    for D in (4//1, 1603//500)
        s = stem(D)
        append!(names, ["results/certificate_$s.json", "results/enclosed_$(s)_L6.json",
            "results/enclosed_$(s)_L12.skeleton.json", "results/enclosed_$(s)_L12.levels.bin"])
    end
    return file_hashes(ROOT, names)
end

"""
    record(row, spec, id)

Serialise a successful seed row; every proof-relevant scalar is an exact rational string.
"""
function record(row, spec, id)
    row.result.bound !== nothing || error("cannot certify a failed seed")
    gap = finite_gap(spec, row.lo, row.hi)
    inf(gap) > 0 || error("finite-ring gap not separated")
    return Dict("id" => id, "reference_D" => string(spec.D),
        "D_interval" => string.([row.lo, row.hi]), "n0" => string(row.n0),
        "beta0" => string(row.beta), "eps_s" => string(row.eps_s),
        "eps_t" => string(row.eps_t), "spatial_lifts" => string(row.lifts),
        "logZ" => endpoints.(row.logZ), "bound" => exact(row.result.bound),
        "finite_gap_L12" => endpoints(gap),
        "ratio_to_gap_L12" => endpoints(box(row.result.bound) / gap),
        "validity" => VALIDITY)
end

"""
    check_record(data, a, b; keep = $KEEP)

Replay one stored case from its spectral inputs. Every stored value must equal the recomputed
one exactly, so any stale or altered field fails, not only an overstated one.
"""
function check_record(data, a, b; keep = KEEP)
    data["validity"] == VALIDITY || error("unsupported validity claim")
    Dref = parse(Q, data["reference_D"])
    Dref == a.D == b.D || error("wrong reference D")
    length(data["D_interval"]) == 2 || error("wrong parameter interval length")
    lo, hi = parse.(Q, data["D_interval"])
    1 <= lo <= hi <= 4 || error("parameter interval outside search range")
    n, lifts = parse(Int, data["n0"]), parse(Int, data["spatial_lifts"])
    0 <= lifts <= 5 && n == 6 * 2^lifts || error("wrong seed length")
    beta = parse(Q, data["beta0"])
    row = evaluate(a, b, beta; lo, hi, max_lifts = lifts, keep)[end]
    length(data["logZ"]) == 4 || error("expected four partition enclosures")
    for (z, stored) in zip(row.logZ, data["logZ"])
        length(stored) == 2 || error("wrong partition interval length")
        zl, zh = parse.(Q, stored)
        (zl, zh) == (Q(inf(z)), Q(sup(z))) || error("log Z enclosure mismatch")
    end
    es, et = parse(Q, data["eps_s"]), parse(Q, data["eps_t"])
    (es, et) == (row.eps_s, row.eps_t) || error("seed defect mismatch")
    result = certified_gap_bound(es, et, beta)
    result.bound !== nothing || error("seed not admissible")
    bound = parse(Q, data["bound"])
    bound == Q(result.bound) > 0 || error("uniform gap mismatch")
    gap = finite_gap(b, lo, hi)
    gl, gh = parse.(Q, data["finite_gap_L12"])
    (gl, gh) == (Q(inf(gap)), Q(sup(gap))) && gl > 0 || error("finite gap mismatch")
    rl, rh = parse.(Q, data["ratio_to_gap_L12"])
    ratio = box(bound) / gap
    (rl, rh) == (Q(inf(ratio)), Q(sup(ratio))) || error("gap ratio mismatch")
    return row
end

"""
    bridge_certificate()

Fresh complete spectral enclosures at D = 4 for the rings L = 8, 10 that the n₀ = 12 seed
does not cover, with their gap intervals rounded outward to a 10⁻⁶ grid.
"""
function bridge_certificate()
    started = time()
    rows = Any[]
    for L in (8, 10)
        # compute_spectrum prints one line per block; the record keeps the block count.
        s = redirect_stdout(() -> compute_spectrum(L, 4//1; rigorous = true), devnull)
        gap = finite_gap((; E = s["energies"], R = s["radii"], L, D = 4//1), 4//1, 4//1)
        inf(gap) > 0 || error("finite gap not separated")
        # The coarse grid absorbs harmless LAPACK/BLAS differences when the checker replays.
        grid = big(10)^6
        lower = floor(BigInt, Q(inf(gap)) * grid) // grid
        higher = ceil(BigInt, Q(sup(gap)) * grid) // grid
        push!(rows, Dict("L" => string(L), "D" => "4//1",
            "gap_interval" => string.([lower, higher]), "dimension" => string(3^L),
            "blocks" => string(length(s["sectors"])), "seconds" => s["seconds"]))
    end
    hashes = proof_hashes()
    return Dict("rings" => rows, "wall_seconds" => time()-started,
        "source_hashes" => hashes,
        "run_metadata" => merge(run_metadata(Dict("D" => "4//1");
            source_sha256 = hashes, packages = (CertifiedGaps, IntervalArithmetic),
            status = "fresh complete spectra; checker repeats these small eigensolves",
            certified = true), Dict("core_sha256" => core_hashes())))
end

"""
    check_bridge(data, threshold)

Recompute the L = 8, 10 enclosures, require each to lie inside its stored interval, and require
the stored lower ends (and the stored L = 6 gap) to be at least `threshold`.
"""
function check_bridge(data, threshold)
    Dict(data["source_hashes"]) == proof_hashes() || error("bridge source mismatch")
    fresh = bridge_certificate()
    length(data["rings"]) == 2 || error("incomplete finite bridge")
    for (stored, recomputed) in zip(data["rings"], fresh["rings"])
        stored["L"] == recomputed["L"] || error("wrong bridge length")
        stored["D"] == recomputed["D"] || error("wrong bridge parameter")
        stored["dimension"] == recomputed["dimension"] || error("wrong bridge dimension")
        stored["blocks"] == recomputed["blocks"] || error("wrong bridge block count")
        sl, sh = parse.(Q, stored["gap_interval"])
        rl, rh = parse.(Q, recomputed["gap_interval"])
        sl <= rl <= rh <= sh || error("bridge did not replay")
        sl >= threshold || error("bridge gap too small")
    end
    a = reference_spectrum(4//1, 6)
    Q(inf(finite_gap(a, 4//1, 4//1))) >= threshold || error("L=6 gap too small")
    return true
end

"""
    check_certificate(path = joinpath(ROOT, "results", "certificate.json"))

Replay every case of the stage-2 certificate and the finite bridge; error on any mismatch.
"""
function check_certificate(path = joinpath(ROOT, "results", "certificate.json"))
    data = parse_exact_json(read(path, String))
    data["format"] == FORMAT || error("format")
    Dict(data["source_hashes"]) == proof_hashes() || error("proof source hash mismatch")
    Dict(data["input_hashes"]) == input_hashes() || error("spectral input hash mismatch")
    parse(Int, data["partition_keep"]) == KEEP || error("unsupported partition cutoff")
    spectra = Dict(D => (reference_spectrum(D, 6), reference_spectrum(D, 12))
        for D in (4//1, 1603//500))
    # Require exactly the advertised cases, so a truncated certificate cannot pass.
    ids = sort([case["id"] for case in data["cases"]])
    ids == sort(EXPECTED_CASES) || error("certificate cases $ids differ from $EXPECTED_CASES")
    for case in data["cases"]
        a, b = spectra[parse(Q, case["reference_D"])]
        row = check_record(case, a, b)
        println(case["id"], ": D in ", case["D_interval"], ", even L >= ", row.n0,
            ", gap >= ", Float64(parse(Q, case["bound"]), RoundDown))
    end
    lifted = only(filter(c -> c["id"] == "space-lift-1-D4//1", data["cases"]))
    bridge = parse_exact_json(read(joinpath(ROOT, "results", "finite_bridge.json"), String))
    check_bridge(bridge, parse(Q, lifted["bound"])) || error("finite bridge failed")
    println("D=4: the n0=12 bound also holds for every even L >= 6 sites (finite bridge)")
    return data
end

end
