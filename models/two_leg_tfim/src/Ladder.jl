"""
    Ladder

The two-leg transverse-field Ising ladder as a core `RingModel` of rungs (d = 4), in the
Hadamard-rotated frame where it is reflection positive and conserves an additive Z₂ charge,
so that the shared core (`GapCertificates.certify_gap`, `VerifiedSpectra.ring_blocks`) does all
of the block reduction, enclosure, log Z and seed arithmetic. Only the model is defined here.
"""
module Ladder
using LinearAlgebra
using CertifiedGaps.RingModels: RingModel, rp_certificate
using CertifiedGaps.VerifiedSpectra: ring_blocks
using CertifiedGaps.Spins: exact_pauli
export POINTS, PROVENANCE, LEG_SWAP, physical_terms, rotated_terms, rung_hadamard, ladder,
    ladder_rp, ladder_energies

const σ = exact_pauli()
const X, Z, I2 = σ.x, σ.z, σ.id

"""
The model-local sources recorded in each certificate's metadata (`Results.model_metadata`)
and enforced by the checker (`Results.verify_model_metadata`).
"""
const PROVENANCE = (dirs = ("src",),
    files = ("certify.jl", "check.jl"))

"Leg exchange a ↔ b on one rung, basis (↑↑, ↑↓, ↓↑, ↓↓)."
const LEG_SWAP = [1 0 0 0; 0 0 1 0; 0 1 0 0; 0 0 0 1]

"Certified parameter points: transverse field g and seed (n₀ rungs, β₀)."
const POINTS = [(; g = 3 // 1, n0 = 4, β0 = 3 // 2)]

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
the Hadamard gate on both spins. W^{⊗L} maps the physical ladder to the rotated one: W X W = Z
spin by spin, so W o_phys W = o_rot and (W ⊗ W) h_phys (W ⊗ W) = h_rot (checked in the tests).
"""
rung_hadamard() = (H = Rational{BigInt}[1 1; 1 -1]; kron(H, H) / 2)

"""
    ladder(g) -> RingModel

The ladder at field g in the rotated frame. The Z₂ charge of a rung basis state is its number
of down spins mod 2 ([0, 1, 1, 0]); it is conserved because every X flips one spin and every
bond or rung term flips two. The leg exchange `LEG_SWAP` is an exact on-site Z₂ symmetry, so
the core also splits every charge and momentum block by its parity.
"""
ladder(g) = RingModel(4, rotated_terms(g).bond, rotated_terms(g).onsite;
    charge = [0, 1, 1, 0], modulus = 2, symmetries = [LEG_SWAP],
    name = "two-leg TFIM ladder, g = $g, rotated frame")

"Exact RP certificate −h = C₁ ⊗ C₁ + C₂ ⊗ C₂ with C₁ = X ⊗ 1, C₂ = 1 ⊗ X (real, R = 1)."
ladder_rp(model) = rp_certificate(model; factors = [kron(X, I2), kron(I2, X)])

"""
    ladder_energies(g, L) -> Vector{Float64}

All 4^L eigenvalues of the ring of L rungs from the core's exact block reduction, by plain
floating-point diagonalization (not verified; for scans and diagnostics only).
"""
function ladder_energies(g, L)
    es = Float64[]
    ring_blocks(ladder(g), L) do M, _, mult, _
        ev = eigvals(Hermitian(M))
        for _ in 1:mult
            append!(es, ev)
        end
    end
    return sort!(es)
end
end
