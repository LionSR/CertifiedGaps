# The dimerized spin-1/2 Heisenberg ring as a core `RingModel` in a sublattice-rotated frame,
# its reflected RP certificate, the seed list and the dimer perturbation bound. Spectra, seed
# defects, gap certificates and all-β exclusions come from `CertifiedGaps` (`VerifiedSpectra`,
# `GapCertificates`, `Exclusions`).
module DimerizedHeisenberg

using LinearAlgebra, IntervalArithmetic
using CertifiedGaps.ExactArithmetic: RationalBig
using CertifiedGaps.RingModels: RingModel, rp_certificate, swap_matrix, is_frame_change
using CertifiedGaps.Recurrence: to_interval
using CertifiedGaps.Spins: exact_spin_half

export DELTAS, SEEDS, DIMER_BETA, MODEL_DIR, FRAME, rotated_model, physical_model, rotated_rp,
    frame_identity, point_params, dimer_bounds

const Q = RationalBig

"The dimerization grid."
const DELTAS = Q[1//2, 3//10, 1//4, 1//5, 3//20, 1//10]

"""
The certified seeds (n₀ cells, β₀) per dimerization, those of the Lean points in
`lean/CertifiedGaps/Models/DimerizedHeisenberg.lean`, plus the prior survey's β₀ = 25/4 at
δ = 1/2. The other grid points have no seed: every β₀ is excluded at n₀ = 2, 4.
"""
const SEEDS = Dict{Q,Vector{Tuple{Int,Q}}}(1//2 => [(4, 6//1), (4, 25//4)],
    3//10 => [(4, 13//2)], 1//4 => [(4, 29//4)])

"β₀ of the recorded (failing) dimer perturbation bound at n₀ = 4: the survey β₀, else 8."
DIMER_BETA(δ) = δ == 1 // 2 ? 25 // big(4) : δ == 3 // 10 ? 13 // big(2) : 8 // big(1)

"The exact parameters of the point δ, recorded in its certificates."
point_params(δ) = Dict("delta" => Q(δ))

const X, IY, Z = let S = exact_spin_half(Q); (S.x, S.iy, S.z) end   # IY = iSʸ, real
const ONE = Matrix{Q}(I, 2, 2)

"The swap of the two spins of a cell: the internal reflection F of the RP certificate."
const SWAP = swap_matrix(2)

"The spin flip σˣ ⊗ σˣ of a cell: a symmetry in both frames."
const SPIN_FLIP = kron(2X, 2X)

"""
The frame rotation W = 1 ⊗ J of one cell, J = 2iSʸ = [0 1; −1 0] on its second spin: the
rotated ring Hamiltonian is U'H U with U = W^{⊗L}.
"""
const FRAME = kron(ONE, 2IY)

# S·S = SˣSˣ + SʸSʸ + SᶻSᶻ, with SʸSʸ = −(iSʸ)(iSʸ).
const HEIS = kron(X, X) - kron(IY, IY) + kron(Z, Z)

"""
    physical_model(δ) -> RingModel

The physical ring of two-spin cells: cell term (1 + δ) S_a·S_b and bond (1 − δ) S_b·S_a'
(second spin of a cell, first spin of the next), with the cell magnetization z_a + z_b as
charge. Used for the frame identity and as an independent Float64 cross-check.
"""
physical_model(δ) = RingModel(4, (1 - δ) * kron(ONE, HEIS, ONE), (1 + δ) * HEIS;
    charge = [1, 0, 0, -1], flip = SPIN_FLIP, reflection = SWAP,
    name = "dimerized Heisenberg, delta = $δ, two-spin cells, physical frame")

# In the rotated frame JᵀSˣJ = −Sˣ, JᵀSᶻJ = −Sᶻ and Jᵀ(iSʸ)J = iSʸ.
rotated(J) = J * (-kron(X, X) - kron(IY, IY) - kron(Z, Z))

"""
    rotated_model(δ) -> RingModel

The ring of two-spin cells in the rotated frame: on-site term (1 + δ)(−SˣSˣ + SʸSʸ − SᶻSᶻ), bond
(1 − δ)(…) between the second spin of a cell and the first of the next. The cell charge
z_a − z_b is conserved and reversed by the cell spin flip. The lattice reflection exchanges the
two spins of every cell, which also reverses the charge, so its on-site part is the cell swap
composed with the spin flip: it preserves the charge and makes every momentum block real.
"""
rotated_model(δ) = RingModel(4, kron(ONE, rotated(1 - δ), ONE), rotated(1 + δ);
    charge = [0, 1, -1, 0], flip = SPIN_FLIP, reflection = SWAP * SPIN_FLIP,
    name = "dimerized Heisenberg, delta = $δ, two-spin cells, rotated frame")

"""
    rotated_rp(δ)

The reflected RP certificate −h = (1 − δ) Σ_α C_α ⊗ F conj(C_α) F, C_α = 1 ⊗ P_α with
P_α ∈ {Sˣ, iSʸ, Sᶻ} and F the cell swap: the hypothesis of the Lean `RPRing` model.
"""
rotated_rp(δ) = rp_certificate(rotated_model(δ); factors = [kron(ONE, P) for P in (X, IY, Z)],
    weights = fill(1 - δ, 3), reflection = SWAP)

"""
    frame_identity(δ) -> Bool

Exact check (`is_frame_change`) that W = 1 ⊗ J is orthogonal and maps the physical cell term
and bond onto those of `rotated_model(δ)`. Since U = W^{⊗L} is a product of cell unitaries,
U'H_L U is the rotated ring Hamiltonian for every L, so both rings have the same spectrum.
"""
frame_identity(δ) = is_frame_change(physical_model(δ), rotated_model(δ), FRAME)

"""
    dimer_bounds(δ, n, β)

Exact norm-perturbation bounds about disconnected strong dimers (δ = 1): with
‖V_n‖ ≤ 3n(1 − δ)/4, every log Z error is ≤ β‖V_n‖. Enclosures of the defects they imply.
"""
function dimer_bounds(delta::Q, n::Int, beta::Q)
    v = 3n * (1 - delta) / 4
    x = exp(-to_interval(beta * (1 + delta)))
    single = (1 + 3x^2) / (1 + 3x)^2
    ps_lower = exp(-to_interval(4beta * v))
    pt_lower = single^(2n) * exp(-to_interval(8beta * v))
    return (; s = 1 - ps_lower, t = 1 - pt_lower, norm_small = v)
end

const MODEL_DIR = normpath(joinpath(@__DIR__, ".."))

end
