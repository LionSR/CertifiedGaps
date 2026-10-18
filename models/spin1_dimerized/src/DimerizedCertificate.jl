"""
    DimerizedCertificate

Gap certificates for the dimerized spin-1 Heisenberg ring
`H_L = Σ_j [1 + (-1)^j δ] S_j · S_{j+1}` (2L spins, L two-spin cells) through the shared core
pipeline. Only the model-specific parts live here: the exact `RingModel` of the two-spin cell in
the sublattice-rotated S^z frame with its lattice reflection, the exact frame unitary, the
reflected RP factors, and the certified points with their seeds (n₀, β₀). Block construction,
eigenvalue enclosure, log Z, seed defects, the rational chain, the drivers `certify_point` /
`check_point` and the floating-point scan spectra are core code.
"""
module DimerizedCertificate
using LinearAlgebra, CertifiedGaps, CertifiedGaps.Results, CertifiedGaps.RingModels
using CertifiedGaps.Spins: integer_spin_one
using CertifiedGaps.GapCertificates: certify_point, check_point, lean_points
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

"The Lean model file whose `points` list holds the seed defects of its theorem (δ = 9/10)."
const LEAN_FILE = normpath(joinpath(ROOT, "..", "..", "lean", "CertifiedGaps", "Models",
    "SpinOneDimerized.lean"))

# Spin-1 in the S^z basis (m = 1, 0, −1) with integer XP = √2 S^x, AP = √2 iS^y (real
# antisymmetric), SZ = S^z, so that S·S = XP⊗XP/2 − AP⊗AP/2 + SZ⊗SZ; P is the flip m ↦ −m.
const XP, AP, SZ, P, ONE = let s = integer_spin_one()
    s.xp, s.ap, s.sz, s.flip, s.id
end
"Rotation by π about y, V = exp(−iπS^y): V S^x V' = −S^x, V S^y V' = S^y, V S^z V' = −S^z."
const V = [0 0 1; 0 -1 0; 1 0 0]
"Swap of the two spins of a cell (state (a, b) has index 3(a−1) + b)."
const SWAP = swap_matrix(3)

"J S·S in the physical S^z basis (exact rationals)."
physical_bond(J) = J * (kron(XP, XP) // 2 - kron(AP, AP) // 2 + kron(SZ, SZ))
"J S·S after the π rotation about y of one of the two spins: −J (XP⊗XP/2 + AP⊗AP/2 + SZ⊗SZ)."
rotated_bond(J) = -J * (kron(XP, XP) // 2 + kron(AP, AP) // 2 + kron(SZ, SZ))
"The frame unitary of one cell, 1 ⊗ V (the second spin b of every cell is rotated)."
cell_unitary() = kron(ONE, V)

"""
    physical(δ) -> RingModel

The ring of two-spin cells (a, b) in the physical S^z basis: on-site term (1 − δ) S_a·S_b, bond
(1 + δ) S_b·S_a' between neighboring cells. Used only for the exact frame check
`is_frame_change(physical(δ), dimerized(δ), cell_unitary())`; the bond is not reflection
positive in this basis.
"""
physical(δ) = RingModel(9, kron(ONE, physical_bond(1 + δ), ONE), physical_bond(1 - δ);
    name = "dimerized spin-1 Heisenberg, delta = $(str(δ)), two-spin cells, physical frame")

"""
    dimerized(δ) -> RingModel

The ring of two-spin cells (a, b) in the rotated frame U = ⊗_cells (1 ⊗ V): on-site term
o = rotated_bond(1 − δ) on (a, b), bond h = rotated_bond(1 + δ) on (b, a') of neighboring cells.
The cell charge z_a − z_b (rotated-frame S^z labels) is conserved, and the flip P ⊗ P maps the
charge sector q onto −q. The lattice reflection is the site reversal of the spin chain: it
reverses the cell order and swaps a ↔ b inside each cell, which maps z_a − z_b to z_b − z_a,
so it is composed with the flip P ⊗ P to keep the charge (on-site part SWAP · (P ⊗ P), checked
exactly by `RingModel`). With it the core builds every momentum block as a real matrix.
"""
function dimerized(δ)
    z = [1, 0, -1]
    RingModel(9, kron(ONE, rotated_bond(1 + δ), ONE), rotated_bond(1 - δ);
        charge = [z[a] - z[b] for a in 1:3 for b in 1:3], flip = kron(P, P),
        reflection = SWAP * kron(P, P),
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

"The seed defects `(s, t)` of the Lean theorem at the point δ, or `nothing` if it has none."
function lean_seed(δ)
    for p in lean_points(LEAN_FILE)
        p.params == [δ, POINTS[δ].beta0] && return (s = p.s, t = p.t)
    end
    return nothing
end

# The model and the keyword arguments that `certify_point` and `check_point` share at δ.
function point(δ)
    (; n0, beta0) = POINTS[δ]
    model = dimerized(δ)
    return model, (; n0, β0 = beta0, params = Dict("delta" => δ), model_dir = ROOT,
        rp = rp(model, δ), positive_transfer = isodd(n0) ? POSITIVE_TRANSFER : "")
end

"""
    certify(δ; progress = false) -> Dict

The certificate of the point δ (`GapCertificates.certify_point` at the seed `POINTS[δ]`, with
the frame recorded), written to `certificate_path(δ)`.
"""
function certify(δ; progress = false)
    model, kw = point(δ)
    record = Dict{String,Any}("frame" =>
        "rotated: U = prod over cells of (1 (x) V), V = exp(-i pi S^y) on spin b")
    return certify_point(certificate_path(δ), model; kw..., record, progress)
end

"""
    check(path; full = false) -> Dict

`GapCertificates.check_point` of a stored certificate against its point (δ from the file name):
model-local source hashes, model, RP data, seed, params, the Lean seed defects where the Lean
file has the point, the exact RP identity and rational chain, and with `full = true` freshly
verified spectra of the rings of n₀ and 2n₀ cells. Throws on any mismatch.
"""
function check(path; full = false)
    m = match(r"^certificate_d(\d+)\.json$", basename(path))
    m === nothing && error("$(basename(path)): not a certificate file name")
    δ = parse(Int, m[1]) // 10
    haskey(POINTS, δ) || error("$(basename(path)): δ = $δ is not a certified point")
    model, kw = point(δ)
    return check_point(path, model; kw..., lean = lean_seed(δ), spectra = full)
end

end
