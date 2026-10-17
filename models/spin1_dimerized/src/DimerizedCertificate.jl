"""
    DimerizedCertificate

Gap certificates for the dimerized spin-1 Heisenberg ring
`H_L = Σ_j [1 + (-1)^j δ] S_j · S_{j+1}` (2L spins, L two-spin cells) through the shared core
pipeline. Only the model-specific parts live here: the exact `RingModel` of the two-spin cell in
the sublattice-rotated S^z frame, its exact frame unitary, its reflected RP factors, the
certified points with their seeds (n₀, β₀), and the loop over `GapCertificates.certify_gap` /
`check_certificate`. Block construction, eigenvalue enclosure, log Z, seed defects and the
rational chain are all core code. The floating-point seed scan lives in `scan_spectra.jl`.
"""
module DimerizedCertificate
using LinearAlgebra, CertifiedGaps, CertifiedGaps.Results, CertifiedGaps.RingModels
using CertifiedGaps.GapCertificates: certify_gap, check_certificate
using CertifiedGaps.ExactArithmetic: parse_rational
BLAS.set_num_threads(1)
# Cheap environment guard: the certificates are produced and replayed with gradual underflow.
get_zero_subnormals() && error("flush-to-zero is enabled; the enclosures assume gradual underflow")

"The model directory (`models/spin1_dimerized`); results live in `ROOT/results`."
const ROOT = normpath(joinpath(@__DIR__, ".."))
"Dimerizations δ of the floating-point scan grid."
const DELTAS = [9//10, 4//5, 7//10, 3//5, 1//2]

"""
Certified points: δ => (n₀ cells, β₀), the best seeds of the floating-point scan
(`scripts/scan.jl`). The odd seed n₀ = 3 rests on the positive cell transfer operator
X_β = K K* ≥ 0 (METHOD.md, "Why the odd seed n₀ = 3 cells is valid").
"""
const POINTS = Dict{Rational{Int},NamedTuple{(:n0, :beta0),Tuple{Int,Rational{Int}}}}(
    9//10 => (n0 = 2, beta0 = 7//2),
    7//10 => (n0 = 3, beta0 = 17//4))

"""
Reference passed to the core as `positive_transfer` for the odd seed n₀ = 3: the cell transfer
operator is X_β = K K* ≥ 0 in the Cartesian frame, and the frame change to the rotated S^z frame
is a product of cell unitaries, under which X_β stays unitarily equivalent.
"""
const POSITIVE_TRANSFER = "models/spin1_dimerized/METHOD.md, \"Why the odd seed n0 = 3 cells is " *
    "valid\": X_beta = K K^* >= 0 (Cartesian frame; the rotated frame differs by cell unitaries)"

"""
Seed defects hardcoded in `lean/CertifiedGaps/Models/SpinOneDimerized.lean` (first entries of
the chain lists `s` and `t`); a certificate must stay below them. δ = 7/10 has no Lean point.
"""
const LEAN_SEEDS = Dict(9//10 => (s = 4972//100000, t = 3206//100000))

# Spin-1 in the S^z basis (m = 1, 0, −1): S^x = XP/√2, iS^y = AP/√2 (real antisymmetric),
# S^z = SZ, so S·S = XP⊗XP/2 − AP⊗AP/2 + SZ⊗SZ with integer XP, AP, SZ.
const ONE = Matrix{Int}(I, 3, 3)
const XP = [0 1 0; 1 0 1; 0 1 0]
const AP = [0 1 0; -1 0 1; 0 -1 0]
const SZ = [1 0 0; 0 0 0; 0 0 -1]
"Rotation by π about y, V = exp(−iπS^y): V S^x V' = −S^x, V S^y V' = S^y, V S^z V' = −S^z."
const V = [0 0 1; 0 -1 0; 1 0 0]
"Spin flip m ↦ −m in the S^z basis: P S^x P = S^x, P S^y P = −S^y, P S^z P = −S^z."
const P = [0 0 1; 0 1 0; 1 0 0]
"Swap of the two spins of a cell (state (a, b) has index 3(a−1) + b)."
const SWAP = let S = zeros(Int, 9, 9)
    for a in 1:3, b in 1:3
        S[3(b - 1) + a, 3(a - 1) + b] = 1
    end
    S
end

"J S·S in the physical S^z basis (exact rationals)."
physical_bond(J) = J * (kron(XP, XP) // 2 - kron(AP, AP) // 2 + kron(SZ, SZ))
"J S·S after the π rotation about y of one of the two spins: −J (XP⊗XP/2 + AP⊗AP/2 + SZ⊗SZ)."
rotated_bond(J) = -J * (kron(XP, XP) // 2 + kron(AP, AP) // 2 + kron(SZ, SZ))
"The frame unitary of one cell, 1 ⊗ V (the second spin b of every cell is rotated)."
cell_unitary() = kron(ONE, V)

"""
    dimerized(δ) -> RingModel

The ring of two-spin cells (a, b) in the rotated frame U = ⊗_cells (1 ⊗ V): on-site term
o = rotated_bond(1 − δ) on (a, b), bond h = rotated_bond(1 + δ) on (b, a') of neighboring cells.
The cell charge z_a − z_b (rotated-frame S^z labels) is conserved, and the flip P ⊗ P maps the
charge sector q onto −q.
"""
function dimerized(δ)
    z = [1, 0, -1]
    RingModel(9, kron(ONE, rotated_bond(1 + δ), ONE), rotated_bond(1 - δ);
        charge = [z[a] - z[b] for a in 1:3 for b in 1:3], flip = kron(P, P),
        name = "dimerized spin-1 Heisenberg, delta = $(str(δ)), two-spin cells, rotated frame")
end

"""
Exact reflected RP certificate −h = Σ_α w_α C_α ⊗ R C_α R with C_α = 1 ⊗ (XP, AP, SZ),
weights (1 + δ)/2, (1 + δ)/2, 1 + δ and R the cell swap.
"""
rp(model, δ) = rp_certificate(model; factors = [kron(ONE, C) for C in (XP, AP, SZ)],
    weights = [(1 + δ) / 2, (1 + δ) / 2, 1 + δ], reflection = SWAP)

"Path of the core-format certificate of δ = d/10."
certificate_path(δ) = joinpath(ROOT, "results", "certificate_d$(Int(10δ)).json")
"Paths of the stored certificates."
result_paths() = sort(filter(p -> occursin(r"^certificate_d\d+\.json$", basename(p)),
    readdir(joinpath(ROOT, "results"); join = true)))
"Format an exact rational as a `\"num//den\"` string."
str(x::Rational) = string(numerator(x), "//", denominator(x))

"""
Model-local sources that produce and check the certificates, as `dirs` and `files` of
`Results.model_metadata` / `verify_model_metadata` (recorded under `source_sha256`, keyed by
repository path; the core and the shared environment `julia/Project.toml`, `julia/Manifest.toml`
are recorded under `core_sha256` as provenance only).
"""
const SOURCES = (dirs = ("src",), files = ("certify.jl", "check.jl"))

"""
    certify_point(δ; progress = false) -> Dict

`certify_gap` at the point's seed `POINTS[δ]`, with δ, the frame and run metadata attached.
"""
function certify_point(δ; progress = false)
    model, (; n0, beta0) = dimerized(δ), POINTS[δ]
    positive_transfer = isodd(n0) ? POSITIVE_TRANSFER : ""
    cert = certify_gap(model; n0, β0 = beta0, rp = rp(model, δ), positive_transfer, progress)
    cert["delta"] = str(δ)
    cert["frame"] = "rotated: U = prod over cells of (1 (x) V), V = exp(-i pi S^y) on spin b"
    config = Dict("delta" => str(δ), "n0" => n0, "beta0" => str(beta0))
    cert["metadata"] = model_metadata(ROOT, config; SOURCES..., packages = (CertifiedGaps,),
        certified = true)
    cert
end

"""
    check_point(path; full = false) -> Dict

Verify one certificate and return the recomputed values; throw on any mismatch. Enforces the
recorded model-local source hashes (and that their key set is exactly the current model-local
set), requires the stored model, RP data and seed to be exactly those of the point, then runs
the core `check_certificate`: with `full = false` the exact RP identity and the exact rational
chain from the stored seed; with `full = true` also freshly verified spectra of the rings of
n₀ and 2n₀ cells, whose defect enclosures must not exceed the stored εₛ, εₜ.
"""
function check_point(path; full = false)
    cert = read_result(path)
    name = basename(path)
    verify_model_metadata(cert["metadata"], ROOT; SOURCES...)
    δ = parse_rational(cert["delta"])
    haskey(POINTS, δ) || error("$name: δ = $δ is not a certified point")
    name == basename(certificate_path(δ)) || error("$name: misnamed")
    model, expected = decode_model(cert["model"]), dimerized(δ)
    model.bond == expected.bond && model.onsite == expected.onsite &&
        model.charge == expected.charge && model.flip == expected.flip ||
        error("$name: stored model differs from the dimerized ring at δ = $δ")
    stored_rp, expected_rp = decode_rp(cert["reflection_positivity"]), rp(expected, δ)
    stored_rp.factors == expected_rp.factors && stored_rp.weights == expected_rp.weights &&
        stored_rp.reflection == expected_rp.reflection ||
        error("$name: stored RP data differ from the reflected Heisenberg factors")
    seed = cert["seed"]
    parse(Int, seed["n0"]) == POINTS[δ].n0 && parse_rational(seed["beta0"]) == POINTS[δ].beta0 ||
        error("$name: stored seed differs from POINTS")
    reference = isodd(POINTS[δ].n0) ? POSITIVE_TRANSFER : ""
    String(get(seed, "positive_transfer", "")) == reference ||
        error("$name: unexpected positive_transfer reference")
    if haskey(LEAN_SEEDS, δ)
        parse_rational(seed["eps_s"]) <= LEAN_SEEDS[δ].s &&
            parse_rational(seed["eps_t"]) <= LEAN_SEEDS[δ].t ||
            error("$name: seed defects exceed the values hardcoded in the Lean statement")
    end
    check_certificate(cert; spectra = full)
end

include("scan_spectra.jl")

end
