# The staggered-field Heisenberg chain on top of the shared core: the exact sublattice rotation
# to the reflection-positive frame, the one-site model and its RP data (the certified claim),
# the two-site blocking with its conserved U(1) cell charge (used only to split the verified seed
# spectra), the certified parameter points, and the hashed source set.
using CertifiedGaps, LinearAlgebra
using CertifiedGaps.RingModels, CertifiedGaps.GapCertificates, CertifiedGaps.Results
using CertifiedGaps.ExactArithmetic: RationalBig, parse_rational
using CertifiedGaps.Spins: exact_spin_half
# The model directory (this file is in src/); results/ and hashed paths are relative to it.
const ROOT = dirname(@__DIR__)

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
a one-site charge here; the U(1) charge needs two-site cells (`BLOCKING`).
"""
site_model(h) = RingModel(2, rotated_bond(), -h * Z; charge = [0, 1], modulus = 2,
    name = "spin-1/2 Heisenberg, staggered field $h, rotated frame")

"Real RP form of the one-site model: −h = Σ C ⊗ C with C ∈ {X, J, Z}, o = −hZ real."
site_rp(h) = rp_certificate(site_model(h); factors = [X, J, Z])

"""
Two-site blocking for the seed spectra. In cells (a, b) with states ↑↑, ↑↓, ↓↑, ↓↓, the rotated
bond flips both spins of a pair together or neither, so the cell charge q = Sᶻ_a − Sᶻ_b
(0, 1, −1, 0) is conserved; it is minus the physical Sᶻ of the cell, since U flips site a. No
on-site cell operator maps q to −q (spin flip reverses the field and the sublattice exchange is
not on-site), so there is no `flip`.
"""
const BLOCKING = (; cell = 2, charge = [0, 1, -1, 0], modulus = 0)

"Model-local sources hashed into every certificate (`model_metadata`, `verify_model_metadata`)."
const SOURCES = (; dirs = ("src",), files = ("certify.jl", "check.jl"))
