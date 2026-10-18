"""
    Ladder

The two-leg transverse-field Ising ladder as a core `RingModel` of rungs (d = 4), in the
Hadamard-rotated frame where it is reflection positive and conserves an additive Z₂ charge,
so that the shared core (`GapCertificates.certify_point`, `VerifiedSpectra.ring_blocks`) does
all of the block reduction, enclosure, log Z and seed arithmetic. Only the model is defined here.
"""
module Ladder
using LinearAlgebra
using CertifiedGaps.RingModels: RingModel, rp_certificate
using CertifiedGaps.Spins: exact_pauli
export POINTS, LEG_SWAP, physical_terms, rotated_terms, rung_hadamard, physical_ladder, ladder,
    ladder_rp, point_path

const σ = exact_pauli()
const X, Z, I2 = σ.x, σ.z, σ.id

"Leg exchange a ↔ b on one rung, basis (↑↑, ↑↓, ↓↑, ↓↓)."
const LEG_SWAP = [1 0 0 0; 0 0 1 0; 0 1 0 0; 0 0 0 1]

"Certified parameter points: transverse field g and seed (n₀ rungs, β₀)."
const POINTS = [(; g = 3 // 1, n0 = 4, β0 = 3 // 2)]

"Certificate file of the point at field g under `root` (results/certificate_g3_1.json)."
point_path(root, g) =
    joinpath(root, "results", "certificate_g$(numerator(g))_$(denominator(g)).json")

"""
    physical_terms(g) -> (; bond, onsite)

Physical-frame rung terms on (leg a ⊗ leg b): h = −(Z_a Z_a′ + Z_b Z_b′) between rungs and
o = −Z_a Z_b − g (X_a + X_b), as in `lean/CertifiedGaps/Models/TwoLegIsing.lean`.
"""
physical_terms(g) = (; bond = -(kron(Z, I2, Z, I2) + kron(I2, Z, I2, Z)),
    onsite = -kron(Z, Z) - g * (kron(X, I2) + kron(I2, X)))

"""
    rotated_terms(g) -> (; bond, onsite)

Rotated-frame terms (Z ↔ X on every spin): h = −(X_a X_a′ + X_b X_b′), o = −X_a X_b − g (Z_a + Z_b).
"""
rotated_terms(g) = (; bond = -(kron(X, I2, X, I2) + kron(I2, X, I2, X)),
    onsite = -kron(X, X) - g * (kron(Z, I2) + kron(I2, Z)))

"""
    rung_hadamard() -> W

The exact rational orthogonal involution W = (H ⊗ H)/2 with H = [1 1; 1 −1] on one rung, i.e.
the Hadamard gate on both spins. W X W = Z spin by spin, so W^{⊗L} maps `physical_ladder(g)`
onto `ladder(g)` (`RingModels.is_frame_change`, checked in the tests).
"""
rung_hadamard() = (H = Rational{BigInt}[1 1; 1 -1]; kron(H, H) / 2)

"The ladder at field g in the physical frame (no charge: the field is off-diagonal there)."
physical_ladder(g) = RingModel(4, physical_terms(g).bond, physical_terms(g).onsite;
    name = "two-leg TFIM ladder, g = $g, physical frame")

"""
    ladder(g) -> RingModel

The ladder at field g in the rotated frame. The Z₂ charge of a rung basis state is its number
of down spins mod 2 ([0, 1, 1, 0]); it is conserved because every X flips one spin and every
bond or rung term flips two. The leg exchange `LEG_SWAP` is an exact on-site Z₂ symmetry, so
the core also splits every charge and momentum block by its parity. The bond is symmetric
under exchanging its two rungs and the reflection maps each rung to itself, so the lattice
reflection x ↦ −x needs no on-site map (`reflection = I`); as the model is real, the core then
builds every momentum block as a real symmetric matrix.
"""
ladder(g) = RingModel(4, rotated_terms(g).bond, rotated_terms(g).onsite;
    charge = [0, 1, 1, 0], modulus = 2, symmetries = [LEG_SWAP],
    reflection = Matrix{Int}(I, 4, 4), name = "two-leg TFIM ladder, g = $g, rotated frame")

"Exact RP certificate −h = C₁ ⊗ C₁ + C₂ ⊗ C₂ with C₁ = X ⊗ 1, C₂ = 1 ⊗ X (real, R = 1)."
ladder_rp(model) = rp_certificate(model; factors = [kron(X, I2), kron(I2, X)])
end
