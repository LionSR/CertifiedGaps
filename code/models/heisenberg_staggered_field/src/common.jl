# The staggered-field Heisenberg chain on top of the shared core: the exact sublattice rotation
# to the reflection-positive frame, the one-site model and its RP data (the certified claim),
# the two-site blocking with its conserved U(1) cell charge (used only to split the verified seed
# spectra), and the certified parameter points.
using CertifiedGaps, LinearAlgebra
using CertifiedGaps.RingModels, CertifiedGaps.GapCertificates, CertifiedGaps.Results
using CertifiedGaps.ExactArithmetic: RationalBig, parse_rational
using CertifiedGaps.Spins: exact_spin_half
# The model directory (this file is in src/); results/ and hashed paths are relative to it.
const ROOT = dirname(@__DIR__)
"The Lean theorem whose `points` carry the seed defects it assumes (read by `lean_points`)."
const LEAN = normpath(joinpath(ROOT, "..", "..", "lean", "CertifiedGaps", "Models",
    "StaggeredHeisenberg.lean"))

"Certified points: staggered field h => seed inverse temperature β₀, at the seed length `N0`."
const SEEDS = [1 // 1 => 21 // 4, 3 // 4 => 21 // 4, 1 // 2 => 11 // 2, 7 // 20 => 23 // 4,
    1 // 4 => 25 // 4]
const FIELDS = first.(SEEDS)
"Seed length n₀ in sites (lengths count sites: the unit cell of the rotated frame is one spin)."
const N0 = 8

tag(h) = replace(string(h), "//" => "_")
certificate_path(h) = joinpath(ROOT, "results", "certificate_$(tag(h)).json")

# Spin-½ matrices X = Sˣ, J = iSʸ (real antisymmetric), Z = Sᶻ, and the identity.
const X, J, Z, I2 = let S = exact_spin_half(Rational{BigInt}); (S.x, S.iy, S.z, S.id) end
"The π rotation about y, U = 2iSʸ = e^{iπSʸ}, applied to every odd site."
const U = 2J

"Physical bond S·S = X⊗X − J⊗J + Z⊗Z (real: SʸSʸ = −J⊗J)."
physical_bond() = kron(X, X) - kron(J, J) + kron(Z, Z)
"Rotated bond −(X⊗X + J⊗J + Z⊗Z): the physical bond conjugated by U on either endpoint."
rotated_bond() = -(kron(X, X) + kron(J, J) + kron(Z, Z))

"""
    site_model(h)

One-site rotated-frame ring: bond −(X⊗X + J⊗J + Z⊗Z) and on-site term −hZ. Unitarily
equivalent (by U on odd sites) to the physical chain on every even ring. Only the Sᶻ parity is
a one-site charge here; the U(1) charge needs two-site cells (`BLOCKING`). The bond is symmetric
under the swap, so the plain lattice reflection (P = I) is a symmetry and is recorded.
"""
site_model(h) = RingModel(2, rotated_bond(), -h * Z; charge = [0, 1], modulus = 2,
    reflection = I2, name = "spin-1/2 Heisenberg, staggered field $h, rotated frame")

"Real RP form of the one-site model: −h = Σ C ⊗ C with C ∈ {X, J, Z}, o = −hZ real."
site_rp(h) = rp_certificate(site_model(h); factors = [X, J, Z])

"""
Two-site blocking for the seed spectra. In cells (a, b) with states ↑↑, ↑↓, ↓↑, ↓↓, the rotated
bond flips both spins of a pair together or neither, so the cell charge q = Sᶻ_a − Sᶻ_b
(0, 1, −1, 0) is conserved; it is minus the physical Sᶻ of the cell, since U flips site a.

The cell model carries neither a `flip` nor a `reflection`: every on-site cell map that mirrors
the bond (the swap of a and b, up to signs) sends q to −q, and the spin flip that would restore
q reverses the field. Diagonalizing the site model instead, with its reflection P = I and its
Sᶻ parity, gives real blocks of about 2070 states at 16 sites, whose enclosure radii give larger
defects than these (charge, momentum) blocks of at most 1620 states (`METHOD.md`).
"""
const BLOCKING = (; cell = 2, charge = [0, 1, -1, 0], modulus = 0)

"""
    point(h)

Keyword arguments of `certify_point` and `check_point` for the field h: the seed (n₀, β₀), the
exact parameter, the model directory, the RP data and the blocking.
"""
point(h) = (; n0 = N0, β0 = Dict(SEEDS)[h], params = Dict("h" => h), model_dir = ROOT,
    rp = site_rp(h), blocking = BLOCKING)
