# Shared definitions of the free-fermion TFIM certificates: the parameter points and seeds, the
# certificate record, and its check. Included by certify.jl, check.jl and test/runtests.jl.

using CertifiedGaps
using CertifiedGaps.ExactArithmetic: RationalBig, frac, round_up, parse_rational
using CertifiedGaps.Results: write_result, read_result, model_metadata, verify_model_metadata
using CertifiedGaps.Seeds: seed_defects
using CertifiedGaps.ReflectionPositivity
using LinearAlgebra
include(joinpath(@__DIR__, "FreeFermions.jl"))

const ROOT = dirname(@__DIR__)
const FORMAT = "CertifiedGaps TFIM free-fermion ring certificate v1"
const PRECISION = 512
# Significant decimal digits kept when the binary defect enclosures are rounded up.
const DIGITS = 8
# Doubling steps of the companion exact chain `rational_gap_chain`. Each step squares the
# defects, so its rounding grid must resolve ε^(2^steps): `chain_grid` takes 40 digits more.
const CHAIN_STEPS = 6
chain_grid(eps_s, eps_t) =
    big(10)^(ceil(Int, -log10(Float64(max(eps_s, eps_t)))) * 2^CHAIN_STEPS + 40)

"""
Parameter points (g, n₀, β₀). The seed lengths are those of the floating-point survey that
listed these points as open targets; β₀ maximizes the interval bound at that n₀ on a grid of
step 1/4 (1/2 for g = 11/10) with scripts/scan.jl.
"""
const POINTS = [
    (g = 2 // 1, n0 = 128, beta0 = 175 // 4),
    (g = 5 // 4, n0 = 144, beta0 = 63 // 1),
    (g = 11 // 10, n0 = 384, beta0 = 359 // 2),
]

tag(g) = "g$(numerator(g))_$(denominator(g))"
certificate_path(g) = joinpath(ROOT, "results", "certificate_$(tag(g)).json.gz")

const PAULI_X = Rational{BigInt}[0 1; 1 0]
const PAULI_Z = Rational{BigInt}[1 0; 0 -1]

"""
    rp_record(g) -> Dict

Exact reflection-positivity data of H = -Σ ZZ - g Σ X in the σᶻ basis: the bond is
h = -Z ⊗ Z, so -h = C ⊗ C̄ with the single real matrix C = Z, and the on-site term o = -gX is
real symmetric. Both identities are checked in exact rational arithmetic here.
"""
function rp_record(g::Rational)
    h = -kron(PAULI_Z, PAULI_Z)
    C = PAULI_Z
    -h == kron(C, conj(C)) || error("-h ≠ Z ⊗ Z̄")
    o = -RationalBig(g) * PAULI_X
    o == transpose(o) && o == conj(o) || error("on-site term is not real symmetric")
    Dict{String,Any}("basis" => "sigma_z", "bond" => "-Z⊗Z", "C" => [["1", "0"], ["0", "-1"]],
        "onsite" => "-g X", "identity" => "-h = Z ⊗ conj(Z)")
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
    seed_enclosures(g, n0, β0) -> NamedTuple

Enclosures of log Z_{n₀}(β₀), log Z_{2n₀}(β₀) and log Z_{2n₀}(2β₀) from the closed form and the
seed defects, rounded up to `DIGITS` significant digits.
"""
function seed_enclosures(g, n0, β0)
    β = RationalBig(β0)
    l = (logZ(n0, g, β; precision = PRECISION), logZ(2n0, g, β; precision = PRECISION),
        logZ(2n0, g, 2β; precision = PRECISION))
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

"""
    certify_point(g, n0, β0) -> Dict

The certificate record of one parameter point (written by certify.jl).
"""
function certify_point(g, n0, β0)
    g, β0 = RationalBig(g), RationalBig(β0)
    s = seed_enclosures(g, n0, β0)
    ib, chain = bounds(s.eps_s, s.eps_t, n0, β0)
    return Dict{String,Any}(
        "format" => FORMAT,
        "model" => "H_L = -sum_j Z_j Z_{j+1} - g sum_j X_j, periodic ring of L sites",
        "g" => string(g), "h" => "0",
        "lengths" => "every even L >= n0 (sites)",
        "exact_gap" => string(2(g - 1)),
        "reflection_positivity" => rp_record(g),
        "precision_bits" => PRECISION,
        "logZ" => Dict("n0_beta0" => endpoints(s.logZ[1]), "2n0_beta0" => endpoints(s.logZ[2]),
            "2n0_2beta0" => endpoints(s.logZ[3])),
        "seed" => Dict{String,Any}("n0" => n0, "beta0" => string(β0),
            "eps_s" => string(s.eps_s), "eps_t" => string(s.eps_t),
            "bound_exact_binary" => string(RationalBig(ib.bound)),
            "interval_bound_decimal" => ib.bound),
        "chain" => chain,
        "metadata" => model_metadata(ROOT, Dict("precision_bits" => PRECISION,
            "digits" => DIGITS); certified = true),
    )
end

"""
    check_point(data, g, n0, β0; full = false)

Check a stored certificate against the definitions above: model-local source hashes, the
point and seed, the exact RP identity, and a rerun of the interval bound and the exact chain
from the stored defects. With `full`, also recompute the log Z enclosures from the closed form
and require the fresh defects to be at most the stored εₛ, εₜ.
"""
function check_point(data, g, n0, β0; full = false)
    g, β0 = RationalBig(g), RationalBig(β0)
    verify_model_metadata(data["metadata"], ROOT)
    data["format"] == FORMAT || error("unknown format")
    parse_rational(data["g"]) == g && data["h"] == "0" || error("parameter mismatch")
    seed = data["seed"]
    seed["n0"] == n0 && parse_rational(seed["beta0"]) == β0 || error("seed mismatch")
    data["reflection_positivity"] == rp_record(g) || error("RP record mismatch")
    eps_s, eps_t = parse_rational(seed["eps_s"]), parse_rational(seed["eps_t"])
    ib, chain = bounds(eps_s, eps_t, n0, β0)
    parse_rational(seed["bound_exact_binary"]) == RationalBig(ib.bound) ||
        error("interval bound mismatch: $(ib.bound)")
    chain["gap_strict_lower"] == data["chain"]["gap_strict_lower"] || error("chain mismatch")
    RationalBig(ib.bound) < 2(g - 1) || error("bound exceeds the exact gap 2(g - 1)")
    if full
        s = seed_enclosures(g, n0, β0)
        s.eps_s <= eps_s && s.eps_t <= eps_t || error("fresh defects exceed the stored ones")
        for (k, x) in zip(("n0_beta0", "2n0_beta0", "2n0_2beta0"), s.logZ)
            data["logZ"][k] == endpoints(x) || error("log Z enclosure $k differs")
        end
    end
    return ib.bound, parse_rational(chain["gap_strict_lower"])
end
