# Shared definitions of the XY-in-a-field certificates: the parameter points and seeds, the
# certificate record, and its check. Included by certify.jl, check.jl and test/runtests.jl.

using CertifiedGaps
using CertifiedGaps.ExactArithmetic: RationalBig, round_up, parse_rational
using CertifiedGaps.Results: write_result, read_result, model_metadata, verify_model_metadata
using CertifiedGaps.Seeds: seed_defects
using LinearAlgebra
include(joinpath(@__DIR__, "FreeFermions.jl"))

const ROOT = dirname(@__DIR__)
const FORMAT = "CertifiedGaps XY-field free-fermion ring certificate v1"
const PRECISION = 512
# Significant decimal digits kept when the binary defect enclosures are rounded up.
const DIGITS = 8
# Doubling steps of the companion exact chain `rational_gap_chain`. Each step squares the
# defects, so its rounding grid must resolve ε^(2^steps): `chain_grid` takes 40 digits more.
const CHAIN_STEPS = 6
chain_grid(eps_s, eps_t) =
    big(10)^(ceil(Int, -log10(Float64(max(eps_s, eps_t)))) * 2^CHAIN_STEPS + 40)

"""
Parameter points (J_x, J_y, h, n₀, β₀): J_x = −1, h = 1 and three values of J_y, with exact gaps
2(h − |J_x + J_y|) = 1/2, 1, 2. β₀ maximizes the interval bound at n₀ = 128 on a grid of step
1/4 (scripts/scan.jl).
"""
const POINTS = [
    (Jx = -1 // 1, Jy = 1 // 4, h = 1 // 1, n0 = 128, beta0 = 189 // 4),
    (Jx = -1 // 1, Jy = 1 // 2, h = 1 // 1, n0 = 128, beta0 = 157 // 4),
    (Jx = -1 // 1, Jy = 1 // 1, h = 1 // 1, n0 = 128, beta0 = 121 // 4),
]

qtag(x) = (x < 0 ? "m" : "") * "$(abs(numerator(x)))_$(denominator(x))"
tag(p) = "Jx$(qtag(p.Jx))_Jy$(qtag(p.Jy))_h$(qtag(p.h))"
certificate_path(p) = joinpath(ROOT, "results", "certificate_$(tag(p)).json.gz")

const CQ = Complex{Rational{BigInt}}
const PAULI_X = CQ[0 1; 1 0]
const PAULI_Y = CQ[0 -im; im 0]
const PAULI_Z = CQ[1 0; 0 -1]

"""
    rp_record(Jx, Jy, h) -> Dict

Exact reflection-positivity data in the standard basis: the bond h_b = J_x X⊗X + J_y Y⊗Y
satisfies −h_b = (−J_x) X ⊗ X̄ + J_y Y ⊗ Ȳ (as Ȳ = −Y), with nonnegative weights for
J_x ≤ 0 ≤ J_y, and the on-site term o = −hZ is real symmetric. Checked here in exact complex
rational arithmetic.
"""
function rp_record(Jx::Rational, Jy::Rational, h::Rational)
    Jx <= 0 <= Jy || error("this frame needs J_x ≤ 0 ≤ J_y")
    bond = Jx * kron(PAULI_X, PAULI_X) + Jy * kron(PAULI_Y, PAULI_Y)
    w = (-RationalBig(Jx), RationalBig(Jy))
    -bond == w[1] * kron(PAULI_X, conj(PAULI_X)) + w[2] * kron(PAULI_Y, conj(PAULI_Y)) ||
        error("RP identity fails")
    o = -RationalBig(h) * PAULI_Z
    o == transpose(o) && o == conj(o) || error("on-site term is not real symmetric")
    Dict{String,Any}("basis" => "standard (sigma_z)", "bond" => "Jx X⊗X + Jy Y⊗Y",
        "C" => ["X", "Y"], "weights" => collect(string.(w)), "onsite" => "-h Z",
        "identity" => "-h_b = (-Jx) X ⊗ conj(X) + Jy Y ⊗ conj(Y)")
end

"Upper bound on the rational `x > 0` with `DIGITS` significant decimal digits."
function round_up_significant(x::RationalBig)
    x > 0 || error("defect enclosure is not positive: $x")
    e = ceil(Int, -log10(Float64(x))) + DIGITS
    return round_up(x, big(10)^e)
end

# Interval enclosure as strings of its BigFloat endpoints (exact binary numbers).
endpoints(x) = [string(exact(inf(x))), string(exact(sup(x)))]

"""
    seed_enclosures(p) -> NamedTuple

Enclosures of log Z_{n₀}(β₀), log Z_{2n₀}(β₀) and log Z_{2n₀}(2β₀) from the closed form and the
seed defects, rounded up to `DIGITS` significant digits.
"""
function seed_enclosures(p)
    β = RationalBig(p.beta0)
    z(L, b) = logZ(L, p.Jx, p.Jy, p.h, b; precision = PRECISION)
    l = (z(p.n0, β), z(2p.n0, β), z(2p.n0, 2β))
    # Same precision as the enclosures: IntervalArithmetic 1.0.12's BigFloat `expm1` returns a
    # point interval that misses the true value when its argument carries more bits than the
    # global precision, and `seed_defects` intersects with that `expm1` enclosure.
    d = setprecision(() -> seed_defects(l...), BigFloat, PRECISION)
    all(x -> precision(inf(x)) == PRECISION && diam(x) > 0, (d.defect_s, d.defect_t)) ||
        error("unexpected precision or width of the defect enclosures")
    # seed_defects converts the upper ends with RationalBig, which rounds to the global
    # precision (PRECISION here); require that conversion to have been exact.
    d.eps_s == exact(sup(d.defect_s)) && d.eps_t == exact(sup(d.defect_t)) ||
        error("inexact conversion of the defect bounds")
    return (; logZ = l, eps_s = round_up_significant(d.eps_s),
        eps_t = round_up_significant(d.eps_t))
end

"The interval bound and the exact chain from the stored defects."
function bounds(eps_s, eps_t, n0, β0)
    ib = certified_gap_bound(eps_s, eps_t, RationalBig(β0))
    ib.bound === nothing && error("the interval recurrence does not certify this seed")
    chain = rational_gap_chain(eps_s, eps_t, n0, RationalBig(β0); steps = CHAIN_STEPS,
        grid = chain_grid(eps_s, eps_t))
    return ib, chain
end

exact_gap(p) = 2(p.h - abs(p.Jx + p.Jy))

"""
    certify_point(p) -> Dict

The certificate record of one parameter point (written by certify.jl).
"""
function certify_point(p)
    s = seed_enclosures(p)
    ib, chain = bounds(s.eps_s, s.eps_t, p.n0, p.beta0)
    return Dict{String,Any}(
        "format" => FORMAT,
        "model" => "H_L = sum_j (Jx X_j X_{j+1} + Jy Y_j Y_{j+1}) - h sum_j Z_j, ring of L sites",
        "Jx" => string(p.Jx), "Jy" => string(p.Jy), "h" => string(p.h),
        "lengths" => "every even L >= n0 (sites)",
        "exact_gap" => string(exact_gap(p)),
        "reflection_positivity" => rp_record(p.Jx, p.Jy, p.h),
        "precision_bits" => PRECISION,
        "logZ" => Dict("n0_beta0" => endpoints(s.logZ[1]), "2n0_beta0" => endpoints(s.logZ[2]),
            "2n0_2beta0" => endpoints(s.logZ[3])),
        "seed" => Dict{String,Any}("n0" => p.n0, "beta0" => string(RationalBig(p.beta0)),
            "eps_s" => string(s.eps_s), "eps_t" => string(s.eps_t),
            "bound_exact_binary" => string(RationalBig(ib.bound)),
            "interval_bound_decimal" => ib.bound),
        "chain" => chain,
        "metadata" => model_metadata(ROOT, Dict("precision_bits" => PRECISION,
            "digits" => DIGITS, "chain_steps" => CHAIN_STEPS); certified = true),
    )
end

"""
    check_point(data, p; full = false)

Check a stored certificate against the point `p`: model-local source hashes, the parameters and
seed, the exact RP record, and a rerun of the interval bound and the exact chain from the stored
defects. With `full`, also recompute the log Z enclosures from the closed form and require the
fresh defects to be at most the stored εₛ, εₜ.
"""
function check_point(data, p; full = false)
    verify_model_metadata(data["metadata"], ROOT)
    data["format"] == FORMAT || error("unknown format")
    all(parse_rational(data[k]) == getfield(p, Symbol(k)) for k in ("Jx", "Jy", "h")) ||
        error("parameter mismatch")
    seed = data["seed"]
    seed["n0"] == p.n0 && parse_rational(seed["beta0"]) == p.beta0 || error("seed mismatch")
    data["reflection_positivity"] == rp_record(p.Jx, p.Jy, p.h) || error("RP record mismatch")
    eps_s, eps_t = parse_rational(seed["eps_s"]), parse_rational(seed["eps_t"])
    ib, chain = bounds(eps_s, eps_t, p.n0, p.beta0)
    parse_rational(seed["bound_exact_binary"]) == RationalBig(ib.bound) ||
        error("interval bound mismatch: $(ib.bound)")
    chain["gap_strict_lower"] == data["chain"]["gap_strict_lower"] || error("chain mismatch")
    RationalBig(ib.bound) < exact_gap(p) || error("bound exceeds the exact gap")
    if full
        s = seed_enclosures(p)
        s.eps_s <= eps_s && s.eps_t <= eps_t || error("fresh defects exceed the stored ones")
        for (k, x) in zip(("n0_beta0", "2n0_beta0", "2n0_2beta0"), s.logZ)
            data["logZ"][k] == endpoints(x) || error("log Z enclosure $k differs")
        end
    end
    return ib.bound, parse_rational(chain["gap_strict_lower"])
end
