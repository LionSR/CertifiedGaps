"""
    LargeDModel

The periodic spin-1 chain H_L = Σⱼ Sⱼ·Sⱼ₊₁ + D Σⱼ (Sᶻⱼ)² as a core `RingModel` in the Sᶻ
eigenbasis (m = +1, 0, −1), its exact reflection-positivity data, the exact change of frame to
the Cartesian basis of the Lean formalization, and the certified parameter points.

With the integer matrices of `Spins.integer_spin_one` (`xp` = √2 Sˣ, `ap` = √2 iSʸ, `sz` = Sᶻ),
S·S = sz⊗sz + xp⊗xp/2 − ap⊗ap/2. The real symmetric involution R = antidiag(1, −1, 1)
(`RP_REFLECTION`) gives R sz R = −sz, R xp R = −xp and R ap R = ap, hence
−h = sz⊗R sz R + ½ xp⊗R xp R + ½ ap⊗R ap R: the "reflected" RP form with integer factors.
The spin flip m ↦ −m reverses the charge and leaves h and o invariant, and the bond is
symmetric under the swap of its two sites, so the lattice reflection x ↦ −x needs no on-site
map (P = 1); with it every momentum block is a real symmetric matrix.
"""
module LargeDModel

using CertifiedGaps.RingModels: RingModel, rp_certificate
using CertifiedGaps.Spins: integer_spin_one

export SPIN, RP_REFLECTION, ring_model, rp_data, cartesian_frame, N0, POINTS, LEAN_FILE,
    STAGE1_SOURCES, stem, spectrum_file, certificate_file

"Spin 1 in the Sᶻ basis with integer entries (`xp`, `ap`, `sz`, `flip`, `id`)."
const SPIN = integer_spin_one()

"The involution R of the reflected RP form; it reverses Sᶻ and Sˣ and keeps Sʸ."
const RP_REFLECTION = [0 0 1; 0 -1 0; 1 0 0]

"The ring at exact anisotropy `D`: Sᶻ charge, spin flip and the lattice reflection P = 1."
ring_model(D::Rational) = RingModel(3,
    kron(SPIN.sz, SPIN.sz) + (kron(SPIN.xp, SPIN.xp) - kron(SPIN.ap, SPIN.ap)) // 2,
    D * SPIN.sz^2; charge = [1, 0, -1], flip = SPIN.flip, reflection = SPIN.id,
    name = "spin-1 chain S.S + D (Sz)^2 in the Sz basis (+1, 0, -1), D = $D")

"Exact RP certificate −h = Σ w C ⊗ R conj(C) R with C ∈ {sz, xp, ap}, w = (1, ½, ½)."
rp_data(model) = rp_certificate(model; factors = [SPIN.sz, SPIN.xp, SPIN.ap],
    weights = [1, 1 // 2, 1 // 2], reflection = RP_REFLECTION)

"""
    cartesian_frame() -> (; W, norms)

Exact change of frame from the Cartesian basis (|x⟩, |y⟩, |z⟩) of the Lean model, where
(S^α)_{ab} = −i ε_{αab}, to the Sᶻ basis used here: the columns of `W` are the Cartesian
coordinates of |+1⟩, |0⟩, |−1⟩, with squared norms `norms`, so U = W diag(norms)^{−1/2} is
unitary. The √2 in U keeps it out of reach of the exact `RingModels.is_frame_change`, so the
identities are asserted for W in Gaussian-integer arithmetic instead: W†W = diag(norms),
W† Sˣ W = xp, W† iSʸ W = ap, W† Sᶻ W = diag(norms) sz and W† (Sᶻ)² W = diag(norms) sz². As
diag(norms)^{1/2} xp diag(norms)^{1/2} = √2 xp (and likewise for ap), these say U† Sˣ U = xp/√2,
U† iSʸ U = ap/√2, U† Sᶻ U = sz and U† (Sᶻ)² U = sz², so (U⊗U)† (S·S) (U⊗U) is the bond
sz⊗sz + xp⊗xp/2 − ap⊗ap/2 of `ring_model`.
"""
function cartesian_frame()
    ε(a, b, c) = (a - b) * (b - c) * (c - a) ÷ 2         # Levi-Civita on {1, 2, 3}
    S = [Complex{Int}[-im * ε(α, a, b) for a in 1:3, b in 1:3] for α in 1:3]
    W = Complex{Int}[-im 0 im; 1 0 1; 0 im 0]
    norms = [2, 1, 2]
    N = [i == j ? norms[i] : 0 for i in 1:3, j in 1:3]
    @assert W' * W == N
    @assert W' * S[1] * W == SPIN.xp && W' * (im * S[2]) * W == SPIN.ap
    @assert W' * S[3] * W == N * SPIN.sz && W' * S[3]^2 * W == N * SPIN.sz^2
    return (; W, norms)
end

"Seed length (sites) of every stage-1 point."
const N0 = 6

"The stage-1 points: exact D and the seed inverse temperature β₀ (as in the Lean file)."
const POINTS = [(D = 4 // 1, beta0 = 9 // 4), (D = 1603 // 500, beta0 = 381 // 125)]

"The Lean theorem whose seed hypotheses (`lean_points`) the stage-1 certificates discharge."
const LEAN_FILE = joinpath(dirname(@__DIR__), "..", "..", "lean", "CertifiedGaps", "Models",
    "SpinOneLargeD.lean")

"""
`model_metadata` sources of the stage-1 certificates: only their producers, so that the
checkers and stage 2 can change without invalidating the multi-hour L = 12 runs.
"""
const STAGE1_SOURCES = (dirs = (), files = ("src/model.jl", "certify.jl"))

"File stem of the stored results for an exact D, e.g. `D1603_500`."
stem(D) = "D$(numerator(D))_$(denominator(D))"

"Path (relative to this model directory) of the stored verified spectrum of the ring L at D."
spectrum_file(D, L) = "results/spectrum_$(stem(D))_L$L.json.gz"

"Path (relative to this model directory) of the stage-1 core certificate at D."
certificate_file(D) = "results/certificate_$(stem(D)).json"

end
