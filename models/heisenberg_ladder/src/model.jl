# The ladder as a core `RingModel` (one rung per cell, d = 4) in the leg-rotated frame, where it
# is reflection positive with the internal cell reflection F = σʸ ⊗ σʸ and conserves the U(1)
# charge 2Sᶻ_total. The physical model is kept only to check the frame unitary exactly.

# Exact spin-½ matrices; A = iSʸ is real antisymmetric, and S·S = X⊗X − A⊗A + Z⊗Z.
const SPIN = exact_spin_half(Q)
const X, A, Z, I2 = SPIN.x, SPIN.iy, SPIN.z, SPIN.id

"Spin-½ exchange S·S on two spins, written with the real matrices X, A = iSʸ, Z."
exchange() = kron(X, X) - kron(A, A) + kron(Z, Z)

"""
Leg bond between rung j (left factor) and rung j+1 (right factor): S·S on each leg. It is the
same in the physical and in the rotated frame, since both of its leg-2 spins are rotated.
"""
function leg_bond()
    leg1(c) = kron(kron(c, I2), kron(c, I2))
    leg2(c) = kron(kron(I2, c), kron(I2, c))
    return leg1(X) - leg1(A) + leg1(Z) + leg2(X) - leg2(A) + leg2(Z)
end

"The frame unitary of one rung: u = exp(iπSʸ) = 2A on the leg-2 spin (real orthogonal)."
frame_unitary() = kron(I2, 2A)

"Cell reflection F = σʸ ⊗ σʸ = −4 A ⊗ A: a real symmetric involution reversing Sᶻ on both legs."
cell_reflection() = -4 * kron(A, A)

"""
Leg exchange in the rotated frame, U · swap · Uᵀ with U the frame unitary of one rung: a real
signed permutation and involution commuting with h, o, the charge and F.
"""
leg_exchange() = frame_unitary() * swap_matrix(2) * transpose(frame_unitary())

"""
    physical_model(J) -> RingModel

The physical ladder in the standard basis: leg bonds S·S, rung term J S₁·S₂ and the charge
2Sᶻ of the rung. It is not reflection positive in this basis; it only checks the frame.
"""
physical_model(J) = RingModel(4, leg_bond(), J * exchange(); charge = [1, 0, 0, -1],
    name = "two-leg Heisenberg ladder, rung coupling $J, physical frame")

"""
    ladder_model(J) -> RingModel

The ladder in the leg-rotated frame U = ⊗_j (1 ⊗ u), u = exp(iπSʸ) on every leg-2 spin: leg
bonds S·S, rung term o = U(J S₁·S₂)U† = −J(X⊗X + A⊗A + Z⊗Z), charge z₁ − z₂ (the physical
total Sᶻ), the charge flip F = σʸ ⊗ σʸ (sectors ±Sᶻ are diagonalized once, Sᶻ = 0 is split
by F-parity) and the on-site Z₂ symmetry `leg_exchange`, by whose parity every block is split.
The leg bond is symmetric under the exchange of its two rungs, so the lattice reflection
j ↦ −j acts with the identity on each rung (`reflection`): k and −k are diagonalized once, and
every block of this real model is a real symmetric matrix.
"""
ladder_model(J) = RingModel(4, leg_bond(), -J * (kron(X, X) + kron(A, A) + kron(Z, Z));
    charge = [0, 1, -1, 0], flip = cell_reflection(), symmetries = [leg_exchange()],
    reflection = Matrix{Int}(I, 4, 4),
    name = "two-leg Heisenberg ladder, rung coupling $J, leg-2 spins rotated by exp(i pi S^y)")

"""
    ladder_rp(model) -> RPCertificate

Exact reflected RP form −h = Σ_C C ⊗ F C F with F = σʸ ⊗ σʸ and the six real factors
C ∈ {X⊗1, A⊗1, Z⊗1, 1⊗X, 1⊗A, 1⊗Z}: F X F = −X, F A F = A, F Z F = −Z on each leg, so the sum is
Σ_legs (−X⊗X + A⊗A − Z⊗Z) = −S·S. Checked by `RingModels.verify_rp`.
"""
ladder_rp(model) = rp_certificate(model;
    factors = vcat([kron(c, I2) for c in (X, A, Z)], [kron(I2, c) for c in (X, A, Z)]),
    reflection = cell_reflection())

"""
    check_frame(J) -> true

Exact identities relating the physical ladder to `ladder_model(J)`: with U the frame unitary
of one rung, Uᵀ maps the physical bond and rung term onto the rotated ones
(`is_frame_change`), and U (Z⊗1 + 1⊗Z) Uᵀ = Z⊗1 − 1⊗Z, so U^{⊗L} maps H_L onto the rotated
ring and the total Sᶻ onto the rotated charge for every L.
"""
function check_frame(J)
    U = exact_matrix(frame_unitary())
    is_frame_change(physical_model(J), ladder_model(J), transpose(U)) ||
        error("the leg rotation does not map the physical ladder onto the rotated one")
    U * exact_matrix(kron(Z, I2) + kron(I2, Z)) * transpose(U) ==
        exact_matrix(kron(Z, I2) - kron(I2, Z)) || error("charge frame fails")
    return true
end
