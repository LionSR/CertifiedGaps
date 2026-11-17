"""
    Spin1SmallSeed

Computer certificate for the leaderboard challenge *a small seed for spin-1*: the spin-1
Heisenberg seed calculation is certified from twisted thermal traces on rings of at most 12
sites (even lengths, the two real twists \$`I`\$ and \$`P`\$, absolute error \$`10^{-3}`\$), lifted
from \$`\\beta = 49/4`\$ to \$`\\beta_0 = 49/2`\$ at moment \$`k = 30`\$ with the stored length-120
integer MPS trial. The exact rational replay (Design C, `real_even`, of
[`code/models/lean_seeds`](../../lean_seeds/README.md)) bounds the seed defects
\$`(p, q) = (\\varepsilon_{\\mathrm{s}}, \\varepsilon_\\tau)`\$ at
\$`(n_0, \\beta_0) = (60, 49/2)`\$ and the rational doubling chain then certifies
\$`\\Delta_L > 0.039036761`\$ for every even \$`L \\ge 60`\$.

The finite trace inputs are inherited from the frozen, hash-audited spin-1 certificate
(`data/models/spin1/source1.json` and `trials.json`, pinned by its `verification.json`):
their rational intervals were produced by exact integer arithmetic and are re-audited here
against the recorded SHA-256 hashes and the challenge's \$`10^{-3}`\$ tolerance. Every scalar
step (sector inversion, polynomial filters, temperature lift, doubling chain) is recomputed
in `Rational{BigInt}` arithmetic by both `certify` and `check`; no floating-point comparison
accepts an inequality.
"""
module Spin1SmallSeed

using JSON3, SHA
using CertifiedGaps.ExactArithmetic
using CertifiedGaps.Results
using CertifiedGaps.RingModels: RingModel, rp_certificate
using CertifiedGaps.Spins: integer_spin_one

export certify, check, MODEL, REPO, RESULTS, SPIN1_DATA, N0, BETA0, EPS, GAP_DECIMAL_DIGITS

# Paths: MODEL is code/models/spin1_small_seed, REPO the repository root.
const MODEL = dirname(@__DIR__)
const CODE = normpath(joinpath(MODEL, "..", ".."))
const REPO = normpath(joinpath(CODE, ".."))
const RESULTS = joinpath(MODEL, "results")
const SPIN1_DATA = joinpath(REPO, "data", "models", "spin1")

# Model-local source sets recorded in the certificate metadata (AGENTS.md "Source hashes").
const META_DIRS = ("src", "test")
const META_FILES = ("certify.jl", "check.jl")

# The per-site energy shift κ = 700741/500000 of the upstream trace tables; it cancels in both
# purities. The traces are Tr[g exp(-β (H_n + n κ I))] on the ring of n spin-1 sites.
const RATSHIFT = frac(700741, 500000)

# Design C ("real_even") parameters: even lengths only, the two real twists, trace error 10⁻³,
# temperature lift 49/4 → 49/2 (multiplier r = 2), moment k = 30, trial length 120.
const EPS = frac(1, 1000)
const NS = (4, 6, 8, 10, 12)
const TWISTS = (0, 1)
const NMAX = 12
const TRIAL = 120
const K = 30
const R = frac(2)
const BASE_BETA = frac(49, 4)
const BETA0 = BASE_BETA * R
const N0 = TRIAL ÷ 2
const GAP_DECIMAL_DIGITS = 9

# Two-twist sector inversion J = (Z + 3P)/4 (twisted sums), N = (Z − P)/4 (non-singlet), with
# per-sector absolute errors eps·Σ|w|; the source conservatively charged eps to every sector.
const JN_WEIGHTS = ([frac(1, 4), frac(3, 4)], [frac(1, 4), frac(-1, 4)])
const JN_MULT = (1, 3)

# Degree-4 squared filters F(x) = x^4 q(x)² selected by the Float64 search in lean_seeds
# (METHOD.md §1, design C); exact_filter re-proves every exterior ray below.
const FILTER_SHIFT = 4
const FILTER_Q = (
    [frac(4359, 500000), frac(0), frac(-32529, 100000), frac(0), frac(1)],
    [frac(11397, 1000000), frac(0), frac(-48363, 125000), frac(0), frac(1)],
)

# Spin 1 in the Sᶻ basis with integer entries; the real symmetric involution
# R = antidiag(1, −1, 1) gives −h = sz⊗R sz R + ½ xp⊗R xp R + ½ ap⊗R ap R exactly.
const SPIN = integer_spin_one()
const RP_REFLECTION = [0 0 1; 0 -1 0; 1 0 0]

"The periodic spin-1 Heisenberg ring as an exact core `RingModel` (Sᶻ basis)."
function ring_model()
    bond = kron(SPIN.sz, SPIN.sz) + (kron(SPIN.xp, SPIN.xp) - kron(SPIN.ap, SPIN.ap)) // 2
    return RingModel(3, bond; charge = [1, 0, -1], flip = SPIN.flip, reflection = SPIN.id,
        name = "spin-1 Heisenberg ring S.S in the Sz basis (+1, 0, -1)")
end

"Exact reflection-positivity certificate; `rp_certificate` verifies it and throws otherwise."
function rp_data()
    factors = [SPIN.sz, SPIN.xp, SPIN.ap]
    weights = [1, 1 // 2, 1 // 2]
    cert = rp_certificate(ring_model(); factors, weights, reflection = RP_REFLECTION)
    return (; cert, factors, weights)
end

"Format the exact rational `q` as a decimal string with `k` fractional digits (q on the 10⁻ᵏ grid)."
function decimal_string(q::RationalBig, k::Int)
    q = round_down(q, big(10)^k)
    n = numerator(q) * fld(big(10)^k, denominator(q))
    s = lpad(string(n), k + 1, '0')
    return s[1:(end - k)] * "." * s[(end - k + 1):end]
end

# ---------------------------------------------------------------------------------------------
# Frozen spin-1 trace inputs: hash audit against the spin-1 certificate's verification.json
# ---------------------------------------------------------------------------------------------

function source1_rows()
    read_exact_json(joinpath(SPIN1_DATA, "source1.json"))["certificates"]
end

function trial_rows()
    read_exact_json(joinpath(SPIN1_DATA, "trials.json"))["certificates"]
end

"""
    audit_spin1_inputs(; full = false) -> Dict

Verify the frozen spin-1 certificate that supplies the finite inputs: `verification.json`
must report `passed`, the SHA-256 of `source1.json` and `trials.json` must equal the
recorded component hashes, and every audited trace row must match the upstream evidence
exactly and carry a rational interval `[lower, upper]` around its center within both its own
`allowed_error` and the challenge tolerance `EPS = 10⁻³`. With `full`, all 18 source-1 rows
are audited; otherwise only the 10 rows of the design (even n ≤ 12, twists I and P).
"""
function audit_spin1_inputs(; full = false)
    verification = read_exact_json(joinpath(SPIN1_DATA, "verification.json"))
    verification["status"] == "passed" ||
        error("spin-1 verification.json status is $(verification["status"])")
    hashes = Dict{String,String}()
    for name in ("source1", "trials")
        path = joinpath(SPIN1_DATA, "$name.json")
        h = sha256_file(path)
        verification["components"][name]["sha256"] == h ||
            error("$path no longer matches the hash recorded in spin-1 verification.json")
        hashes[name] = h
    end
    rows = source1_rows()
    selected = [r for r in rows if parse(Int, r["n"]) in NS && parse(Int, r["twist"]) in TWISTS]
    length(selected) == length(NS) * length(TWISTS) ||
        error("source1.json does not contain the $(length(NS) * length(TWISTS)) design rows")
    audited = full ? rows : selected
    for row in audited
        row["matches_upstream_exactly"] == true ||
            error("trace row n=$(row["n"]) twist=$(row["twist"]) does not match upstream")
        c = parse_rational(row["center"])
        lo = parse_rational(row["lower"])
        hi = parse_rational(row["upper"])
        allowed = parse_rational(row["allowed_error"])
        c - allowed < lo <= hi < c + allowed ||
            error("trace row n=$(row["n"]) twist=$(row["twist"]) interval violates its error bound")
        c - EPS < lo <= hi < c + EPS ||
            error("trace row n=$(row["n"]) twist=$(row["twist"]) exceeds the 10⁻³ tolerance")
    end
    return Dict("hashes" => hashes, "rows_audited" => length(audited),
        "selected" => [Dict("n" => parse(Int, r["n"]), "twist" => parse(Int, r["twist"]),
            "total" => r["total"], "center" => r["center"], "lower" => r["lower"],
            "upper" => r["upper"], "allowed_error" => r["allowed_error"]) for r in selected])
end

"Twisted-trace centers Z_n(g) as exact rationals, keyed by (n, twist)."
function stored_table()
    Dict((parse(Int, row["n"]), parse(Int, row["twist"])) => parse_rational(row["center"])
        for row in source1_rows())
end

"""
    verify_trial() -> Dict

The stored exact length-120 integer MPS trial: its energy per bond `U/V` beats the shift,
`500000 U + 700741 V < 0` with `V > 0`, hash-audited as in `audit_spin1_inputs`.
"""
function verify_trial()
    row = only(r for r in trial_rows() if parse(Int, r["length"]) == TRIAL)
    u, v = parse(BigInt, row["U"]), parse(BigInt, row["V"])
    v > 0 && 500000 * u + 700741 * v < 0 ||
        error("the stored length-$TRIAL trial does not beat the shift")
    return Dict("length" => TRIAL, "U" => string(u), "V" => string(v),
        "energy_per_bond" => string(u // v), "positive_margin" => string(-500000 * u - 700741 * v),
        "scope" => "stored exact trial output of data/models/spin1/trials.json; " *
            "no fresh contraction")
end

# ---------------------------------------------------------------------------------------------
# Exact scalar chain: sector moments, filters, temperature lift, doubling recurrence
# ---------------------------------------------------------------------------------------------

"Sector trace moments: J_n = (Z_n + 3 P_n)/4 and N_n = (Z_n − P_n)/4 at every length."
function sector_moments(table)
    ns = sort(unique(first(k) for k in keys(table)))
    [Dict(n => sum(w[g + 1] * table[(n, g)] for g in 0:(length(w) - 1)) for n in ns)
        for w in JN_WEIGHTS]
end

const SECTOR_ERRORS = map(w -> EPS * sum(abs, w), JN_WEIGHTS)

function convsquare(q, shift)
    out = zeros(eltype(q), 2length(q) - 1 + shift)
    for i in eachindex(q), j in eachindex(q)
        out[i + j - 1 + shift] += q[i] * q[j]
    end
    out
end

"Exact exterior-ray check of the squared filter F(x) = x^shift q(x)² against mass bound B."
function filter_valid(P, B, v)
    for sign in (-1, 1)
        shifted = [sum(binomial(i, j) * P[i + 1] * sign^i * v^(i - j) for i in j:(length(P) - 1))
            for j in 0:(length(P) - 1)]
        first(shifted) > B && all(>=(0), shifted[2:end]) || return false
    end
    true
end

"Filter modulus cap: the exact least v for which both exterior rays are valid."
function exact_filter(mom, eps, shift, q; digits = 10)
    @assert iseven(shift) && shift >= 4
    P = convsquare(q, shift)
    B = sum(c * mom[i - 1] for (i, c) in enumerate(P) if c != 0) + eps * sum(abs, P)
    @assert B > 0
    _, cap = first_true(v -> filter_valid(P, B, v), frac(1, 20), frac(2); digits)
    @assert filter_valid(P, B, cap)
    return Dict("shift" => string(shift), "q" => string.(q), "P" => string.(P), "B" => string(B),
        "cap" => string(cap))
end

"Capped-mass upper bound: count·cap^r + (m − count·cap)^r with count = ⌊m/cap⌋."
function capped_upper(m, cap, r)
    @assert m >= 0 && cap > 0 && r >= 1
    count = fld(m, cap)
    count * power_upper(cap, r) + power_upper(m - count * cap, r)
end

"Sector mass bounds at moment k, and the combined M = Σ mult·b, E = b_J (JN mode)."
function exact_bounds(masses, filters, k)
    r = frac(k, NMAX)
    bs = RationalBig[]
    for (m, f) in zip(masses, filters)
        push!(bs, capped_upper(m, parse_rational(f["cap"])^NMAX, r))
    end
    sum(JN_MULT .* bs), bs[1]
end

"Residual convex bound at the two vertices of the JN simplex, in exact rationals."
function exact_residual(D, tauj, taui, t, s)
    @assert 0 < t <= taui <= tauj <= D && s >= 1
    vertices = ((frac(0), frac(0)), (tauj - t, frac(0)))
    maximum(power_upper(x, s) + 2 * power_upper(y / 2, s) + 3 * power_upper((D - t - x - y) / 3, s)
        for (x, y) in vertices)
end

"""
    exact_seed(masses, filters) -> Dict

Temperature lift at moment `k = K`, multiplier `r = R`, trial length `TRIAL`: exact upper
bounds `p ≥ 1 − P_s(n0, β0)` and `q ≥ 1 − P_τ(2n0, β0)` at `(n0, β0) = (60, 49/2)`.
"""
function exact_seed(masses, filters)
    k, r = K, R
    @assert TRIAL % 4 == 0 && iseven(k) && NMAX <= k <= TRIAL ÷ 2 && r >= 1
    M, E = exact_bounds(masses, filters, k)
    witness_grid = big(10)^12
    D = round_up(power_upper(M, r), witness_grid)
    E = round_up(power_upper(E, r), witness_grid)
    h = frac(TRIAL, 2k)
    tauj = (D + 3E) / 4
    taui = tauj
    t = round_down(first(power_interval(inv(D), inv(2h - 1))), witness_grid) -
        inv(frac(witness_grid))
    @assert 0 < t < taui && 3t > D
    @assert D * power_upper(t, 2h - 1) < 1
    initial = t
    refinements = Any[]
    for _ in 1:35
        residual = exact_residual(D, tauj, taui, t, 2h)
        @assert residual < 1
        next = round_down(first(power_interval(1 - residual, inv(2h))), witness_grid) -
               inv(frac(witness_grid))
        next <= t && break
        @assert next < taui && power_upper(next, 2h) + residual < 1
        push!(refinements, Dict("old" => string(t), "new" => string(next),
            "residual" => string(residual)))
        t = next
    end
    res = exact_residual(D, tauj, taui, t, h)
    dom = first(power_interval(t, h))
    p = 1 - (1 + res / dom)^(-2)
    Q, _ = exact_bounds(masses, filters, TRIAL)
    Q = round_up(Q, witness_grid)
    @assert Q > 1
    ratio = round_up(power_upper(Q - 1, r), big(10)^18)
    q = 1 - (1 + ratio)^(-2)
    return Dict("n0" => N0, "beta0" => string(BETA0), "k" => k, "r" => string(r), "h" => string(h),
        "trial" => TRIAL, "D" => string(D), "E" => string(E), "tauj" => string(tauj),
        "taui" => string(taui), "initial_t" => string(initial), "refinements" => refinements,
        "t" => string(t), "Q" => string(Q), "spatial_residual" => string(res),
        "dominant_power" => string(dom), "physical_ratio" => string(ratio),
        "p" => string(p), "q" => string(q),
        "diagonal_threshold_passes" => max(p, q) < frac(12731, 100000))
end

zq(a, b) = 1 - (1 - b) * (1 - a)^2
fq(z) = z^2 / (2 * (1 - z)^2)

"Exact replay of the doubling-chain validity predicate (`ChainCertificate.Valid`)."
function chain_valid(c)
    s, t = parse_rational.(c["s"]), parse_rational.(c["t"])
    r, epsilon = parse_rational(c["r"]), parse_rational(c["epsilon"])
    J = length(s) - 1
    length(t) == J + 1 && r > 0 && 0 < epsilon <= frac(1, 100) || return false
    s[end] <= epsilon && t[end] <= epsilon && r^(2^J) * (5epsilon) <= 1 || return false
    all(x -> 0 <= x <= 1, vcat(s, t)) || return false
    for j in 1:J
        z = zq(s[j], t[j]); zt = zq(t[j], s[j + 1])
        z < 1 && zt < 1 && fq(z) <= s[j + 1] && fq(zt) <= t[j + 1] || return false
        r^(2^(j - 1)) * z < 1 - z || return false
    end
    true
end

"""
    exact_chain(p, q, beta; steps = 9) -> Dict

Nine exact doubling steps on the 10⁻⁵¹² rounding grid, the chain rate, and the strict gap
lower bound `Δ_L > log(r)/β` for every even `L ≥ n0`.
"""
function exact_chain(p, q, beta; steps = 9)
    grid = big(10)^512
    s, t = [p], [q]
    for _ in 1:steps
        z = zq(last(s), last(t)); @assert 0 < z < frac(1, 2)
        p = round_up(fq(z), grid)
        zt = zq(last(t), p); @assert 0 < zt < 1
        q = round_up(fq(zt), grid)
        push!(s, p); push!(t, q)
    end
    epsilon = max(last(s), last(t)); @assert epsilon <= frac(1, 100)
    validrate(r) = r^(2^steps) * (5epsilon) <= 1 &&
        all(r^(2^(j - 1)) * zq(s[j], t[j]) < 1 - zq(s[j], t[j]) for j in 1:steps)
    @assert validrate(frac(1)) && !validrate(frac(10))
    r, _ = first_true(x -> !validrate(x), frac(1), frac(10); digits = 8)
    @assert r > 1
    c = Dict("s" => string.(s), "t" => string.(t), "r" => string(r), "epsilon" => string(epsilon))
    @assert chain_valid(c)
    lower = first(log_enclosure(r)) / beta
    claim = round_down(lower, big(10)^9) - frac(1, big(10)^9)
    @assert 0 < claim < lower
    c["beta"] = string(beta)
    c["gap_strict_lower"] = string(claim)
    c["gap_display"] = Float64(claim)
    return c
end

# ---------------------------------------------------------------------------------------------
# Certificate writer and checker
# ---------------------------------------------------------------------------------------------

"The certificate's exact core, computed from scratch; shared by `certify` and `check`."
function compute()
    audit = audit_spin1_inputs()
    table = stored_table()
    moments = sector_moments(table)
    filters = [exact_filter(m, e, FILTER_SHIFT, q)
        for (m, e, q) in zip(moments, SECTOR_ERRORS, FILTER_Q)]
    masses = [m[NMAX] + e for (m, e) in zip(moments, SECTOR_ERRORS)]
    trial = verify_trial()
    seed = exact_seed(masses, filters)
    p, q = parse_rational(seed["p"]), parse_rational(seed["q"])
    chain = exact_chain(p, q, BETA0)
    # The chain dominates the original published target log(20)/784 without floats.
    ratio = frac(784) / BETA0
    @assert parse_rational(chain["r"])^numerator(ratio) > frac(20)^denominator(ratio)
    rp = rp_data()
    return Dict("audit" => audit,
        "moments" => [Dict(string(n) => string(v) for (n, v) in m) for m in moments],
        "filters" => filters, "masses" => string.(masses), "trial" => trial, "seed" => seed,
        "chain" => chain,
        "rp" => Dict("factors" => [[string(x) for x in f] for f in rp.factors],
            "weights" => string.(rp.weights),
            "identity" => "-h = sz⊗R sz R + 1/2 xp⊗R xp R + 1/2 ap⊗R ap R, " *
                "verified exactly by rp_certificate"))
end

"Assemble and write `results/certificate.json`."
function certify()
    core = compute()
    seed, chain = core["seed"], core["chain"]
    claim = chain["gap_strict_lower"]
    gap_decimal = decimal_string(parse_rational(claim), GAP_DECIMAL_DIGITS)
    cert = Dict{String,Any}(
        "format" => "CertifiedGaps ring gap certificate v1 (spin-1 small seed, rings <= 12 sites)",
        "model" => Dict("name" => "spin-1 Heisenberg ring",
            "hamiltonian" => "H_L = sum_j S_j . S_{j+1} on a periodic ring of L sites",
            "shift" => string(RATSHIFT),
            "note" => "traces are Tr[g exp(-beta (H_n + n*shift))] with g in {I, P}"),
        "seed" => Dict("n0" => N0, "beta0" => string(BETA0), "eps_s" => seed["p"],
            "eps_t" => seed["q"], "k" => seed["k"], "r" => seed["r"], "trial" => seed["trial"],
            "gap_strict_lower" => claim),
        "interval_bound_decimal" => gap_decimal,
        "gap_lower_bound" => claim,
        "traces" => Dict("source" => "data/models/spin1/source1.json",
            "sha256" => core["audit"]["hashes"]["source1"], "eps" => string(EPS),
            "rows" => core["audit"]["selected"]),
        "trial" => merge(core["trial"],
            Dict("sha256" => core["audit"]["hashes"]["trials"])),
        "input" => Dict("mode" => "JN", "nmax" => NMAX, "eps" => string(EPS),
            "base_beta" => string(BASE_BETA), "masses" => core["masses"],
            "filters" => core["filters"]),
        "chain" => chain,
        "rp" => core["rp"],
        "lengths" => Dict("L0" => N0, "parity" => "even", "boundary" => "periodic"),
        "metadata" => model_metadata(MODEL,
            "spin-1 small seed: even rings <= 12 sites, real twists, trace error 10^-3, " *
            "lift 49/4 -> 49/2, moment 30, trial 120"; dirs = META_DIRS, files = META_FILES,
            packages = (JSON3,), certified = true),
    )
    path = joinpath(RESULTS, "certificate.json")
    write_result(path, cert)
    println("wrote $path")
    println("seed (n0, beta0) = ($N0, $BETA0), defects p = $(Float64(parse_rational(seed["p"]))), ",
        "q = $(Float64(parse_rational(seed["q"]))), chain rate r = ",
        "$(Float64(parse_rational(chain["r"]))), gap > $gap_decimal for even L >= $N0")
    return path
end

"""
    check(; full = false)

Recompute the certificate from the stored inputs and compare every field with
`results/certificate.json`; throws (fails loudly) on any mismatch. With `full`, all 18
source-1 trace rows are re-audited instead of only the 10 design rows.
"""
function check(; full = false)
    path = joinpath(RESULTS, "certificate.json")
    stored = read_exact_json(path)
    verify_model_metadata(stored["metadata"], MODEL; dirs = META_DIRS, files = META_FILES)
    core = compute()
    full && audit_spin1_inputs(; full = true)
    # Exact RP identity: rp_certificate throws unless -h = Σ w C ⊗ R conj(C) R holds exactly.
    rp_data()
    # read_exact_json returns bare integers as strings; compare parsed values.
    parse(Int, stored["seed"]["n0"]) == N0 && stored["seed"]["beta0"] == string(BETA0) ||
        error("stored seed does not match the design (n0, beta0) = ($N0, $BETA0)")
    stored["seed"]["eps_s"] == core["seed"]["p"] && stored["seed"]["eps_t"] == core["seed"]["q"] ||
        error("stored seed defects differ from the recomputed ones")
    stored["seed"]["gap_strict_lower"] == core["chain"]["gap_strict_lower"] ||
        error("stored gap bound differs from the recomputed chain bound")
    stored["interval_bound_decimal"] ==
        decimal_string(parse_rational(core["chain"]["gap_strict_lower"]), GAP_DECIMAL_DIGITS) ||
        error("stored decimal bound is not the rounded-down exact bound")
    stored["chain"] == core["chain"] || error("stored chain differs from the replay")
    stored["input"]["masses"] == core["masses"] || error("stored masses differ")
    stored["input"]["filters"] == core["filters"] || error("stored filters differ")
    stored["trial"]["U"] == core["trial"]["U"] && stored["trial"]["V"] == core["trial"]["V"] ||
        error("stored trial differs")
    rows = stored["traces"]["rows"]
    length(rows) == length(NS) * length(TWISTS) || error("stored trace rows are incomplete")
    table = stored_table()
    for row in rows
        key = (parse(Int, row["n"]), parse(Int, row["twist"]))
        haskey(table, key) || error("stored trace row $key is not in source1.json")
        table[key] == parse_rational(row["center"]) || error("stored trace center for $key differs")
    end
    stored["traces"]["sha256"] == core["audit"]["hashes"]["source1"] ||
        error("stored trace hash differs")
    # The validator-facing fields: the certificate states seed (n0, beta0) with the bound.
    @assert chain_valid(stored["chain"])
    p = parse_rational(stored["seed"]["eps_s"]); q = parse_rational(stored["seed"]["eps_t"])
    @assert max(p, q) < frac(12731, 100000)
    claim = parse_rational(stored["chain"]["gap_strict_lower"])
    @assert 0 < claim < first(log_enclosure(parse_rational(stored["chain"]["r"]))) / BETA0
    println("certificate checked", full ? " (full input audit)" : "",
        ": seed ($N0, $BETA0), gap > $(stored["interval_bound_decimal"]) for even L >= $N0")
    return true
end

end
