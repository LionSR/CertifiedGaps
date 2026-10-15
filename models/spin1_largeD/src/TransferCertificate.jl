module TransferCertificate

include("../model.jl")
include("../spectra.jl")

const ROOT = dirname(@__DIR__)
const Q = Rational{BigInt}
const box = CertifiedGaps.to_interval
const VALIDITY = "unique ground state and gap >= bound for every even L >= n0 sites"
exact(x) = string(Q(x))
endpoints(x) = [exact(inf(x)), exact(sup(x))]
upper(x) = box(max(0, sup(x)))

"Load all verified levels, checking the historical certificate's byte fingerprint."
function reference_spectrum(D, L)
    stem = "D$(numerator(D))_$(denominator(D))"
    cert = parse_exact_json(read(joinpath(ROOT, "results", "certificate_$stem.json"), String))
    name = "enclosed_$(stem)_L$L"
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

"Enclose log Z using verified low levels and an explicit positive omitted-level tail."
function logz(spec, beta; keep = 8192)
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
        floor = box(E[k+1]) - box(maximum(@view R[k+1:end]))
        tail = (length(E) - k) * exp(-b * (floor - shift))
        total += interval(0.0, sup(tail))
    end
    return log(total) - b * shift
end

"Centered Weyl transport: log Z for H(D)-L(D-Dref)/2 on an exact D interval."
function transported_logz(z, L, beta, Dref, lo, hi)
    lo <= hi || throw(ArgumentError("reversed D interval"))
    radius = sup(box(L) * box(beta) * box(max(abs(lo-Dref), abs(hi-Dref))) / 2)
    return z + interval(-radius, radius)
end

"The purity-squaring bound, capped by the trivial defect bound one."
function square_defect(s)
    0 <= inf(s) <= sup(s) <= 1 || throw(ArgumentError("invalid defect"))
    sup(s) == 1 && return interval(1.0)
    return box(min(1, sup(s^2 / (2 * (1-s)^2))))
end

"Initial and space-lifted seeds from four log-partition enclosures at L=6 and L=12."
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

"Build every requested seed from a reference spectrum and an exact parameter interval."
function evaluate(a, b, beta; lo = a.D, hi = lo, max_lifts = 5, keep = 8192)
    a.L == 6 && b.L == 12 && a.D == b.D || error("expected matching L=6,12 spectra")
    rp_identity(lo) && rp_identity(hi) || error("RP identity failed")
    base = [logz(a, beta; keep), logz(b, beta; keep),
        logz(a, 2beta; keep), logz(b, 2beta; keep)]
    return evaluate_logs(base, beta, a.D, lo, hi; max_lifts)
end

"Reuse finite log Z enclosures while scanning parameter intervals or seed lengths."
function evaluate_logs(base, beta, Dref, lo, hi; max_lifts = 5)
    zs = [transported_logz(z, L, temp, Dref, lo, hi) for (z, L, temp) in
        zip(base, (6, 12, 6, 12), (beta, beta, 2beta, 2beta))]
    rows = seed_ladder(zs...; max_lifts)
    failed = (; bound = nothing, upper = nothing, j_admissible = nothing, log_rJ = nothing)
    return [(; row..., beta, lo, hi, logZ = zs,
        result = row.eps_s < 1 && row.eps_t < 1 &&
            (1-row.eps_s)^2 * (1-row.eps_t) > 1//2 ?
            certified_gap_bound(row.eps_s, row.eps_t, beta) : failed) for row in rows]
end

"Enclose the L=12 gap at every D in the parameter interval, including Weyl error."
function finite_gap(spec, lo, hi)
    width = box(spec.L) * box(max(abs(lo-spec.D), abs(hi-spec.D)))
    # Blockwise enclosures can cross after sorting the combined list. A common radius
    # encloses the globally sorted eigenvalues without assuming their block identities.
    r = 2 * box(maximum(spec.R)) + width
    return box(spec.E[2]) - box(spec.E[1]) + interval(-sup(r), sup(r))
end

"Files that define the proof arithmetic and its environment, including shared-core sources."
function proof_hashes()
    hashes = source_hashes(ROOT; dirs = ("src",),
        files = ("model.jl", "spectra.jl", "transfer_certify.jl", "transfer_check.jl",
            "bridge.jl", "test/transfer.jl", "Project.toml", "Manifest.toml"))
    core = source_hashes(joinpath(ROOT, "../../julia"); dirs = ("src",),
        files = ("Project.toml",))
    merge!(hashes, Dict("../../julia/$f" => h for (f, h) in core))
    return hashes
end

"Byte hashes of the original complete spectral evidence and its certificate anchors."
function input_hashes()
    names = String[]
    for D in (4//1, 1603//500)
        stem = "D$(numerator(D))_$(denominator(D))"
        append!(names, ["results/certificate_$stem.json", "results/enclosed_$(stem)_L6.json",
            "results/enclosed_$(stem)_L12.skeleton.json",
            "results/enclosed_$(stem)_L12.levels.bin"])
    end
    return file_hashes(ROOT, names)
end

"Serialize a successful interval seed; every proof-relevant scalar is an exact string."
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

"Replay one case from its spectral inputs; reject changed seeds, scope, or bounds."
function check_record(data, a, b; keep = 8192)
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
        zl <= Q(inf(z)) <= Q(sup(z)) <= zh || error("log Z enclosure mismatch")
    end
    es, et = parse(Q, data["eps_s"]), parse(Q, data["eps_t"])
    row.eps_s <= es < 1 && row.eps_t <= et < 1 || error("invalid seed upper bounds")
    result = certified_gap_bound(es, et, beta)
    result.bound !== nothing || error("seed not admissible")
    bound = parse(Q, data["bound"])
    0 < bound <= Q(result.bound) || error("overstated uniform gap")
    gap = finite_gap(b, lo, hi)
    gl, gh = parse.(Q, data["finite_gap_L12"])
    0 < gl <= Q(inf(gap)) <= Q(sup(gap)) <= gh || error("finite gap mismatch")
    rl, rh = parse.(Q, data["ratio_to_gap_L12"])
    ratio = box(bound) / gap
    rl <= Q(inf(ratio)) <= Q(sup(ratio)) <= rh || error("gap ratio mismatch")
    return row
end

end
