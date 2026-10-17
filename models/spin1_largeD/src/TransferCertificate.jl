"""
    TransferCertificate

Stage 2 of the spin-1 large-D certificate, built on the core stage-1 outputs: the verified
spectra of the rings L = 6 and L = 12 that `certify.jl` stores with
`CertifiedGaps.VerifiedSpectra.write_spectrum` (SHA-256 pinned in the stage-1 certificates).

From log Z enclosures at a reference D (`log_partition_enclosure`), this module lifts the seed
length by transfer-operator squaring (n₀ = 6·2ʲ) with the exact purity coupling identity,
transports seeds to exact D intervals by a centered Weyl bound, and passes the exact rational
defects to `CertifiedGaps.certified_gap_bound`. It also builds and replays the finite L = 8, 10
bridge at D = 4. The inequalities are proved in `METHOD.md`.
"""
module TransferCertificate

using CertifiedGaps, IntervalArithmetic
using CertifiedGaps.Results: write_result, read_exact_json, model_metadata, verify_model_metadata,
    sha256_file
using CertifiedGaps.VerifiedSpectra: verified_ring_spectrum, log_partition_enclosure,
    spectral_gap_enclosure, read_spectrum
include("model.jl")
using .LargeDModel

const ROOT = dirname(@__DIR__)
const Q = Rational{BigInt}
const box = CertifiedGaps.to_interval
const FORMAT = "spin1_largeD transfer and Weyl certificate v2"
const VALIDITY = "unique ground state and gap >= bound for every even L >= n0 sites"
"The case ids `scripts/transfer_certify.jl` writes; `check_certificate` requires exactly these."
const EXPECTED_CASES = ["reference-D4//1", "survey-D4-beta3",
    ["space-lift-$k-D4//1" for k in 1:4]..., "transport-D4//1",
    "reference-D1603//500", "transport-D1603//500"]
"""
Lower ends of the transported D intervals, fixed across regenerations so that the certified
claims stay comparable; `scripts/transfer_certify.jl` reports by bisection how far the current
enclosures would reach (a diagnostic only).
"""
const TRANSPORT_LO = Dict(4 // 1 => 39979773 // 10000000, 1603 // 500 => 16029999 // 5000000)
exact(x) = string(Q(x))
endpoints(x) = [exact(inf(x)), exact(sup(x))]
upper(x) = box(max(0, sup(x)))

"""
    reference_spectrum(D, L)

The stored verified spectrum of the ring of `L` sites at `D`, after checking its bytes against
the SHA-256 recorded in the stage-1 certificate and its size against 3^L.
"""
function reference_spectrum(D, L)
    cert = read_exact_json(joinpath(ROOT, certificate_file(D)))
    file = spectrum_file(D, L)
    hash = sha256_file(joinpath(ROOT, file))
    hash == cert["spectrum_sha256"][file] || error("spectrum hash mismatch for $file")
    spectrum = read_spectrum(joinpath(ROOT, file))
    sum(length(s.centers) * s.multiplicity for s in spectrum) == 3^L ||
        error("wrong spectral dimension in $file")
    return (; D, L, spectrum, hash)
end

"Enclosure of log Z_L(β) from a stored verified spectrum."
logz(spec, beta) = log_partition_enclosure(spec.spectrum, beta)

"""
    transported_logz(z, L, beta, Dref, lo, hi)

Centered Weyl transport: enclosure of log Z for H(D) − L(D−Dref)/2 at every D in [lo, hi],
from an enclosure `z` at Dref. The constant shift cancels in both purities.
"""
function transported_logz(z, L, beta, Dref, lo, hi)
    lo <= hi || throw(ArgumentError("reversed D interval"))
    radius = sup(box(L) * box(beta) * box(max(abs(lo - Dref), abs(hi - Dref))) / 2)
    return z + interval(-radius, radius)
end

"The purity-squaring bound s ↦ s²/(2(1−s)²), capped by the trivial defect bound one."
function square_defect(s)
    0 <= inf(s) <= sup(s) <= 1 || throw(ArgumentError("invalid defect"))
    sup(s) == 1 && return interval(1.0)
    return box(min(1, sup(s^2 / (2 * (1 - s)^2))))
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
        t = upper(1 - (1 - s2) * (1 - t)^2)
        push!(rows, (; n0 = 6 * 2^j, eps_s = Q(sup(s)), eps_t = Q(sup(t)), lifts = j))
    end
    return rows
end

"""
    evaluate(a, b, beta; lo = a.D, hi = lo, max_lifts = 5)

All seeds of the ladder for reference spectra `a` (L = 6) and `b` (L = 12), transported to
the exact interval [lo, hi]. See `evaluate_logs`.
"""
function evaluate(a, b, beta; lo = a.D, hi = lo, max_lifts = 5)
    a.L == 6 && b.L == 12 && a.D == b.D || error("expected matching L=6,12 spectra")
    base = [logz(a, beta), logz(b, beta), logz(a, 2beta), logz(b, 2beta)]
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
                      (1 - row.eps_s)^2 * (1 - row.eps_t) > 1 // 2
    return [(; row..., beta, lo, hi, logZ = zs, result = admissible(row) ?
        certified_gap_bound(row.eps_s, row.eps_t, beta) : failed) for row in rows]
end

"""
    finite_gap(spec, lo, hi) -> (lower, upper)

Exact bounds on the gap of the ring at every D in [lo, hi]: the core's `spectral_gap_enclosure`
widened by the Weyl shift L·max|D − D_ref| (each of the two lowest levels of the centered
Hamiltonian moves by at most half of it), rounded outward to binary64 values so that the stored
rationals stay short. `lower > 0` proves a unique ground state.
"""
function finite_gap(spec, lo, hi)
    lower, higher = spectral_gap_enclosure(spec.spectrum)
    width = spec.L * max(abs(lo - spec.D), abs(hi - spec.D))
    return Q(inf(box(Q(lower) - width))), Q(sup(box(Q(higher) + width)))
end

"Float64 interval enclosing an exact (lower, upper) pair."
gap_box((lo, hi)) = interval(inf(box(lo)), sup(box(hi)))

"""
    record(row, spec, id)

Serialize a successful seed row; every proof-relevant scalar is an exact rational string.
"""
function record(row, spec, id)
    row.result.bound !== nothing || error("cannot certify a failed seed")
    gap = finite_gap(spec, row.lo, row.hi)
    gap[1] > 0 || error("finite-ring gap not separated")
    return Dict("id" => id, "reference_D" => string(spec.D),
        "D_interval" => string.([row.lo, row.hi]), "n0" => string(row.n0),
        "beta0" => string(row.beta), "eps_s" => string(row.eps_s),
        "eps_t" => string(row.eps_t), "spatial_lifts" => string(row.lifts),
        "logZ" => endpoints.(row.logZ), "bound" => exact(row.result.bound),
        "finite_gap_L12" => string.(collect(gap)),
        "ratio_to_gap_L12" => endpoints(box(Q(row.result.bound)) / gap_box(gap)),
        "validity" => VALIDITY)
end

"""
    check_record(data, a, b)

Replay one stored case from its spectral inputs. Every stored value must equal the recomputed
one exactly, so any stale or altered field fails, not only an overstated one.
"""
function check_record(data, a, b)
    data["validity"] == VALIDITY || error("unsupported validity claim")
    Dref = parse(Q, data["reference_D"])
    Dref == a.D == b.D || error("wrong reference D")
    length(data["D_interval"]) == 2 || error("wrong parameter interval length")
    lo, hi = parse.(Q, data["D_interval"])
    1 <= lo <= hi <= 4 || error("parameter interval outside search range")
    n, lifts = parse(Int, data["n0"]), parse(Int, data["spatial_lifts"])
    0 <= lifts <= 5 && n == 6 * 2^lifts || error("wrong seed length")
    beta = parse(Q, data["beta0"])
    row = evaluate(a, b, beta; lo, hi, max_lifts = lifts)[end]
    length(data["logZ"]) == 4 || error("expected four partition enclosures")
    for (z, stored) in zip(row.logZ, data["logZ"])
        length(stored) == 2 || error("wrong partition interval length")
        Tuple(parse.(Q, stored)) == (Q(inf(z)), Q(sup(z))) || error("log Z enclosure mismatch")
    end
    es, et = parse(Q, data["eps_s"]), parse(Q, data["eps_t"])
    (es, et) == (row.eps_s, row.eps_t) || error("seed defect mismatch")
    result = certified_gap_bound(es, et, beta)
    result.bound !== nothing || error("seed not admissible")
    bound = parse(Q, data["bound"])
    bound == Q(result.bound) > 0 || error("uniform gap mismatch")
    gap = finite_gap(b, lo, hi)
    Tuple(parse.(Q, data["finite_gap_L12"])) == gap && gap[1] > 0 ||
        error("finite gap mismatch")
    data["ratio_to_gap_L12"] == endpoints(box(bound) / gap_box(gap)) ||
        error("gap ratio mismatch")
    return row
end

"Model-local sources of the stage-2 proof, hashed by `metadata` and enforced by the checker."
const STAGE2_SOURCES = (dirs = ("src",),
    files = ("scripts/transfer_certify.jl", "check.jl"))

"`model_metadata` of a stage-2 result (model-local `STAGE2_SOURCES` and core hashes)."
metadata(config) = model_metadata(ROOT, config; STAGE2_SOURCES...,
    packages = (CertifiedGaps, IntervalArithmetic), certified = true)

"SHA-256 of the stage-1 certificates and spectra that stage 2 reads."
function input_hashes()
    names = [f for (D, _) in POINTS for f in
             (certificate_file(D), spectrum_file(D, 6), spectrum_file(D, 12))]
    return Dict(name => sha256_file(joinpath(ROOT, name)) for name in names)
end

"""
    bridge_certificate()

Fresh verified spectra at D = 4 for the rings L = 8, 10 that the n₀ = 12 seed does not cover,
with the core gap lower bounds rounded down to a 10⁻⁶ grid.
"""
function bridge_certificate()
    started = time()
    rows = Any[]
    for L in (8, 10)
        spectrum = verified_ring_spectrum(ring_model(4 // 1), L)
        gap = Q(spectral_gap_enclosure(spectrum)[1])
        gap > 0 || error("finite gap not separated")
        # The coarse grid absorbs harmless LAPACK/BLAS differences when the checker replays.
        lower = floor(BigInt, gap * big(10)^6) // big(10)^6
        push!(rows, Dict("L" => string(L), "D" => "4//1", "gap_lower" => string(lower),
            "dimension" => string(3^L), "blocks" => string(length(spectrum)),
            "seconds" => sum(s.label["seconds"] for s in spectrum)))
    end
    return Dict("rings" => rows, "wall_seconds" => time() - started,
        "metadata" => metadata(Dict("D" => "4//1")))
end

"""
    check_bridge(data, threshold)

Recompute the L = 8, 10 spectra, require each fresh gap lower bound to be at least the stored
one, and require the stored ones (and the stored L = 6 gap) to be at least `threshold`.
"""
function check_bridge(data, threshold)
    verify_model_metadata(data["metadata"], ROOT; STAGE2_SOURCES...)
    fresh = bridge_certificate()
    length(data["rings"]) == 2 || error("incomplete finite bridge")
    for (stored, recomputed) in zip(data["rings"], fresh["rings"])
        for key in ("L", "D", "dimension", "blocks")
            stored[key] == recomputed[key] || error("bridge $key mismatch")
        end
        sl, rl = parse(Q, stored["gap_lower"]), parse(Q, recomputed["gap_lower"])
        sl <= rl || error("bridge did not replay")
        sl >= threshold || error("bridge gap too small")
    end
    finite_gap(reference_spectrum(4 // 1, 6), 4 // 1, 4 // 1)[1] >= threshold ||
        error("L=6 gap too small")
    return true
end

"""
    check_certificate(path = joinpath(ROOT, "results", "certificate.json"))

Replay every case of the stage-2 certificate and the finite bridge; error on any mismatch.
"""
function check_certificate(path = joinpath(ROOT, "results", "certificate.json"))
    data = read_exact_json(path)
    data["format"] == FORMAT || error("format")
    verify_model_metadata(data["metadata"], ROOT; STAGE2_SOURCES...)
    Dict(data["input_hashes"]) == input_hashes() || error("spectral input hash mismatch")
    spectra = Dict(D => (reference_spectrum(D, 6), reference_spectrum(D, 12))
                   for (D, _) in POINTS)
    # Require exactly the advertised cases, so a truncated certificate cannot pass.
    ids = sort([case["id"] for case in data["cases"]])
    ids == sort(EXPECTED_CASES) || error("certificate cases $ids differ from $EXPECTED_CASES")
    for case in data["cases"]
        a, b = spectra[parse(Q, case["reference_D"])]
        row = check_record(case, a, b)
        startswith(case["id"], "transport-") &&
            parse(Q, case["D_interval"][1]) != TRANSPORT_LO[a.D] &&
            error("transported interval differs from the advertised one")
        println(case["id"], ": D in ", case["D_interval"], ", even L >= ", row.n0,
            ", gap >= ", Float64(parse(Q, case["bound"]), RoundDown))
    end
    lifted = only(filter(c -> c["id"] == "space-lift-1-D4//1", data["cases"]))
    bridge = read_exact_json(joinpath(ROOT, "results", "finite_bridge.json"))
    check_bridge(bridge, parse(Q, lifted["bound"])) || error("finite bridge failed")
    println("D=4: the n0=12 bound also holds for every even L >= 6 sites (finite bridge)")
    return data
end

end
