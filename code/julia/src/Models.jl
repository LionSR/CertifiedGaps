"""
    CertifiedGaps.Models

Ring models used to test the criterion. Each constructor returns a named tuple
`(; chain, rp_bond, d)`: the `Chain` used for exact diagonalisation, and the bond term in the
local basis where reflection positivity holds (check it with `ReflectionPositivity.rp_test`).
"""
module Models

using LinearAlgebra
using ..Spins: spin_matrices, cartesian_generators, PAULI
using ..ExactDiagonalization: Chain, ring_bonds

export tfim, heisenberg, staggered_heisenberg, tfim_logZ, tfim_seed_defects

"Transverse-field Ising ring −Σ ZZ − g Σ X (RP in the σᶻ basis: −h = Z ⊗ Z̄)."
function tfim(g)
    h = -kron(PAULI.z, PAULI.z)
    return (chain = Chain(2, ring_bonds(h); onsite = -g * PAULI.x), rp_bond = h, d = 2)
end

"""
    heisenberg(s)

Antiferromagnetic spin-s Heisenberg ring Σ Sⱼ·Sⱼ₊₁, with Sᶻ sectors (charge 2m). `rp_bond` is
the same bond in the Cartesian basis, −S·S = Σ_α G_α ⊗ G_α, available for s = 1, 2; for s = 1
this is the model of the note.
"""
function heisenberg(s)
    S = spin_matrices(s)
    bond = real(kron(S.sz, S.sz) + (kron(S.sp, S.sm) + kron(S.sm, S.sp)) / 2)
    charge = [Int(2 * (s - k)) for k in 0:Int(2s)]
    rp_bond = s in (1, 2) ? -sum(kron(G, G) for G in cartesian_generators(Int(s))) : nothing
    return (chain = Chain(Int(2s + 1), ring_bonds(bond); charge), rp_bond, d = Int(2s + 1))
end

"""
    staggered_heisenberg(hs)

Spin-½ Heisenberg ring with a staggered field, Σ Sⱼ·Sⱼ₊₁ − h_s Σ (−1)ʲ Sᶻⱼ (even L). ED uses
the physical frame with Sᶻ sectors, the field folded into the bond to the right of each site.
After the sublattice rotation e^{iπSʸ} on odd sites every bond is −SˣSˣ + SʸSʸ − SᶻSᶻ (RP with
real C = {Sˣ, iSʸ, Sᶻ}) and the field becomes the uniform real on-site term −h_s Sᶻ; for even L
the two rings are unitarily equivalent.
"""
function staggered_heisenberg(hs)
    S = spin_matrices(1 // 2)
    heis = real(kron(S.sz, S.sz) + (kron(S.sp, S.sm) + kron(S.sm, S.sp)) / 2)
    field = kron(S.sz, I(2))
    bonds(L) = [(j, mod1(j + 1, L), heis - hs * (-1)^j * field) for j in 1:L]
    rotated = real(-kron(S.sx, S.sx) + kron(S.sy, S.sy) - kron(S.sz, S.sz))
    return (chain = Chain(2, bonds; charge = [1, -1]), rp_bond = rotated, d = 2)
end

# Exact partition function of the transverse-field Ising ring, g > 1. Jordan–Wigner with the
# parity projection (survey/models/tfim_g2/run.py, §4): with ε_k = 2√(1 + g² − 2g cos k),
# antiperiodic momenta A = {2π(m + ½)/L} and periodic P = {2πm/L},
#     Z_L(β) = ½ [Π_A 2cosh(βε_k/2) + Π_A 2sinh(βε_k/2) + Π_P 2cosh(βε_k/2) − Π_P 2sinh(βε_k/2)].
# The gap tends to 2(g − 1).

_tfim_energies(L, g, shift) = [2sqrt(1 + g^2 - 2g * cospi(2(m + shift) / L)) for m in 0:(L-1)]

"""
    tfim_logZ(L, β, g)

log Z_L(β) of the TFIM ring in `BigFloat` (current precision), in the cancellation-free form
log Z = −log 2 + log A⁺ + log(1 + A⁻/A⁺ + (P⁺/A⁺)(1 − P⁻/P⁺)).
"""
function tfim_logZ(L::Integer, β, g)
    β, g = BigFloat(β), BigFloat(g)
    function logs(ε)                        # (log Π 2cosh(βε/2), log Π 2sinh(βε/2))
        x = β .* ε ./ 2
        return sum(x .+ log1p.(exp.(-2x))), sum(x .+ log1p.(-exp.(-2x)))
    end
    Ap, Am = logs(_tfim_energies(L, g, big(1) / 2))
    Pp, Pm = logs(_tfim_energies(L, g, big(0)))
    return -log(big(2)) + Ap + log(1 + exp(Am - Ap) - exp(Pp - Ap) * expm1(Pm - Pp))
end

"Exact seed defects (1 − Pₛ(n, β), 1 − P_τ(2n, β)) of the TFIM ring, as `BigFloat`s."
function tfim_seed_defects(n::Integer, β, g; precision = 2048)
    setprecision(BigFloat, precision) do
        lZ(L, b) = tfim_logZ(L, b, g)
        -expm1(lZ(2n, β) - 2lZ(n, β)), -expm1(lZ(2n, 2β) - 2lZ(2n, β))
    end
end

end
