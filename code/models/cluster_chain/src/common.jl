# Shared definitions of the cluster-chain certificates: the parameter point and seed, the exact
# reflection-positivity data in the rotated cell frame, the certificate record, and its check.
# Included by certify.jl, check.jl and test/runtests.jl.

using CertifiedGaps
using CertifiedGaps.ExactArithmetic: RationalBig, round_up, parse_rational
using CertifiedGaps.Results: write_result, read_result, model_metadata, verify_model_metadata
using CertifiedGaps.Seeds: seed_defects
using CertifiedGaps.RingModels: RingModel, rp_certificate, verify_rp, is_frame_change
using LinearAlgebra
include(joinpath(@__DIR__, "FreeFermions.jl"))

const ROOT = dirname(@__DIR__)
const FORMAT = "CertifiedGaps cluster-chain free-fermion ring certificate v1"
const PRECISION = 512
# Significant decimal digits kept when the binary defect enclosures are rounded up.
const DIGITS = 8
# Doubling steps of the companion exact chain `rational_gap_chain`. Each step squares the
# defects, so its rounding grid must resolve ε^(2^steps): `chain_grid` takes 40 digits more.
const CHAIN_STEPS = 6
chain_grid(eps_s, eps_t) =
    big(10)^(ceil(Int, -log10(Float64(max(eps_s, eps_t)))) * 2^CHAIN_STEPS + 40)

"""
Parameter point (g, n₀, β₀), with n₀ in two-site cells: the open target g = 1/2 at the survey's
n₀ = 128 cells; β₀ maximizes the interval bound on a grid of step 1/4 (scripts/scan.jl).
"""
const POINTS = [(g = 1 // 2, n0 = 128, beta0 = 175 // 2)]

tag(g) = "g$(numerator(g))_$(denominator(g))"
certificate_path(g) = joinpath(ROOT, "results", "certificate_$(tag(g)).json.gz")

const CQ = Complex{RationalBig}
const PX = CQ[0 1; 1 0]
const PZ = CQ[1 0; 0 -1]
const I2 = CQ[1 0; 0 1]

"""
    physical_model(g) -> RingModel

The cluster chain on two-site cells (2c, 2c+1): the bond between neighboring cells carries
Z_{2c+1} X_{2c+2} Z_{2c+3} and Z_{2c} X_{2c+1} Z_{2c+2}, the cell term −g(X ⊗ 1 + 1 ⊗ X).
"""
function physical_model(g)
    bond = -(kron(kron(I2, PZ), kron(PX, PZ)) + kron(kron(PZ, PX), kron(PZ, I2)))
    return RingModel(4, bond, -RationalBig(g) * (kron(PX, I2) + kron(I2, PX)); name = "cluster")
end

"The cell frame U = (H ⊗ H) CS, with H ⊗ H = ½[±1] exact and CS = diag(1, 1, 1, i)."
function cell_frame()
    H2 = CQ[1 1; 1 -1]
    return kron(H2, H2) / 2 * CQ[1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 im]
end

"The model in the cell frame: bond (U ⊗ U)† h (U ⊗ U) and cell term U† o U."
function rotated_model(g)
    m, U = physical_model(g), cell_frame()
    UU = kron(U, U)
    return RingModel(4, UU' * m.bond * UU, U' * m.onsite * U; name = "cluster, rotated")
end

entry_string(z) = imag(z) == 0 ? string(real(z)) : "$(real(z)) + $(imag(z))im"
matrix_strings(A) = [[entry_string(A[i, j]) for j in axes(A, 2)] for i in axes(A, 1)]

"""
    rp_record(g) -> Dict

Exact RP data: the frame change from the physical to the rotated model (core
`is_frame_change`), the factorization −h' = Σ w_α C_α ⊗ C̄_α of the rotated bond found by the
core's exact LDL† (`rp_certificate`), and its exact check (`verify_rp`, which also requires the
rotated cell term to be real).
"""
function rp_record(g)
    phys, rot, U = physical_model(g), rotated_model(g), cell_frame()
    is_frame_change(phys, rot, U) || error("U is not an exact frame change")
    cert = rp_certificate(rot)
    verify_rp(rot, cert) || error("RP certificate fails")
    return Dict{String,Any}("frame" => "U = (H⊗H) CS, CS = diag(1, 1, 1, i)",
        "U" => matrix_strings(U), "weights" => string.(cert.weights),
        "factors" => [matrix_strings(C) for C in cert.factors],
        "identity" => "-h' = sum_a w_a C_a ⊗ conj(C_a), h' = (U⊗U)' h (U⊗U)")
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

Enclosures of log Z_{n₀}(β₀), log Z_{2n₀}(β₀) and log Z_{2n₀}(2β₀) (n₀ in cells) from the
closed form and the seed defects, rounded up to `DIGITS` significant digits.
"""
function seed_enclosures(g, n0, β0)
    β = RationalBig(β0)
    z(L, b) = logZ(L, g, b; precision = PRECISION)
    l = (z(n0, β), z(2n0, β), z(2n0, 2β))
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
    g = RationalBig(g)
    s = seed_enclosures(g, n0, β0)
    ib, chain = bounds(s.eps_s, s.eps_t, n0, β0)
    return Dict{String,Any}(
        "format" => FORMAT,
        "model" => "H = -sum_j Z_{j-1} X_j Z_{j+1} - g sum_j X_j, ring of L two-site cells",
        "g" => string(g),
        "lengths" => "every even L >= n0 (two-site cells)",
        "exact_gap" => string(2abs(1 - g)),
        "reflection_positivity" => rp_record(g),
        "precision_bits" => PRECISION,
        "logZ" => Dict("n0_beta0" => endpoints(s.logZ[1]), "2n0_beta0" => endpoints(s.logZ[2]),
            "2n0_2beta0" => endpoints(s.logZ[3])),
        "seed" => Dict{String,Any}("n0_cells" => n0, "beta0" => string(RationalBig(β0)),
            "eps_s" => string(s.eps_s), "eps_t" => string(s.eps_t),
            "bound_exact_binary" => string(RationalBig(ib.bound)),
            "interval_bound_decimal" => ib.bound),
        "chain" => chain,
        "metadata" => model_metadata(ROOT, Dict("precision_bits" => PRECISION,
            "digits" => DIGITS, "chain_steps" => CHAIN_STEPS); certified = true),
    )
end

"""
    check_point(data, g, n0, β0; full = false)

Check a stored certificate against the point: model-local source hashes, the parameter and
seed, the exact RP record (recomputed), and a rerun of the interval bound and the exact chain
from the stored defects. With `full`, also recompute the log Z enclosures from the closed form
and require the fresh defects to be at most the stored εₛ, εₜ.
"""
function check_point(data, g, n0, β0; full = false)
    g, β0 = RationalBig(g), RationalBig(β0)
    verify_model_metadata(data["metadata"], ROOT)
    data["format"] == FORMAT || error("unknown format")
    parse_rational(data["g"]) == g || error("parameter mismatch")
    seed = data["seed"]
    seed["n0_cells"] == n0 && parse_rational(seed["beta0"]) == β0 || error("seed mismatch")
    data["reflection_positivity"] == rp_record(g) || error("RP record mismatch")
    eps_s, eps_t = parse_rational(seed["eps_s"]), parse_rational(seed["eps_t"])
    ib, chain = bounds(eps_s, eps_t, n0, β0)
    parse_rational(seed["bound_exact_binary"]) == RationalBig(ib.bound) ||
        error("interval bound mismatch: $(ib.bound)")
    chain["gap_strict_lower"] == data["chain"]["gap_strict_lower"] || error("chain mismatch")
    RationalBig(ib.bound) < 2abs(1 - g) || error("bound exceeds the exact gap")
    if full
        s = seed_enclosures(g, n0, β0)
        s.eps_s <= eps_s && s.eps_t <= eps_t || error("fresh defects exceed the stored ones")
        for (k, x) in zip(("n0_beta0", "2n0_beta0", "2n0_2beta0"), s.logZ)
            data["logZ"][k] == endpoints(x) || error("log Z enclosure $k differs")
        end
    end
    return ib.bound, parse_rational(chain["gap_strict_lower"])
end
