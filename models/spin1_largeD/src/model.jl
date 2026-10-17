"""
    LargeDModel

The periodic spin-1 chain H_L = Σⱼ Sⱼ·Sⱼ₊₁ + D Σⱼ (Sᶻⱼ)² as a core `RingModel` in the Sᶻ
eigenbasis (m = +1, 0, −1), its exact reflection-positivity data, the exact unitary to the
Cartesian frame of the Lean formalization, and the certified parameter points.

In this basis S⁺ = √2 `SP`, so S·S = Sᶻ⊗Sᶻ + SP⊗SM + SM⊗SP has integer entries. The reflection
R = antidiag(1, −1, 1) maps Sᶻ ↦ −Sᶻ and S^± ↦ −S^∓, hence
−h = Sᶻ⊗R Sᶻ R + SP⊗R SP R + SM⊗R SM R: the "reflected" RP form with integer factors. The spin
flip m ↦ −m (a permutation) reverses the charge and leaves h and o invariant.
"""
module LargeDModel

using CertifiedGaps.RingModels: RingModel, rp_certificate

export SZ, SP, SM, REFLECTION, FLIP, ring_model, rp_data, cartesian_frame, N0, POINTS,
    LEAN_SEED_BOUNDS, STAGE1_SOURCES, stem, spectrum_file, certificate_file

const SZ = [1 0 0; 0 0 0; 0 0 -1]
const SP = [0 1 0; 0 0 1; 0 0 0]          # S⁺/√2
const SM = permutedims(SP)                 # S⁻/√2
const REFLECTION = [0 0 1; 0 -1 0; 1 0 0]
const FLIP = [0 0 1; 0 1 0; 1 0 0]

"The ring at exact anisotropy `D`, with the Sᶻ charge and the spin flip."
ring_model(D::Rational) = RingModel(3, kron(SZ, SZ) + kron(SP, SM) + kron(SM, SP),
    D * [1 0 0; 0 0 0; 0 0 1]; charge = [1, 0, -1], flip = FLIP,
    name = "spin-1 chain S.S + D (Sz)^2 in the Sz basis (+1, 0, -1), D = $D")

"Exact RP certificate −h = Σ C ⊗ R conj(C) R with C ∈ {Sᶻ, SP, SM} and R = `REFLECTION`."
rp_data(model) = rp_certificate(model; factors = [SZ, SP, SM], reflection = REFLECTION)

"""
    cartesian_frame() -> (; W, norms)

Exact change of frame from the Cartesian basis (|x⟩, |y⟩, |z⟩) of the Lean model, where
(S^α)_{ab} = −i ε_{αab}, to the Sᶻ basis used here: the columns of `W` are the Cartesian
coordinates of |+1⟩, |0⟩, |−1⟩, with squared norms `norms`, so U = W diag(norms)^{−1/2} is
unitary. Asserts in exact Gaussian-integer arithmetic that W†W = diag(norms),
W† Sᶻ W = diag(norms) SZ, W† S^± W = 2 SP, 2 SM and W† (Sᶻ)² W = diag(norms) diag(1, 0, 1),
i.e. U† Sᶻ U = SZ, U† S^± U = √2 SP, √2 SM, U†(Sᶻ)²U = diag(1, 0, 1). Then
(U⊗U)† (S·S) (U⊗U) = SZ⊗SZ + SP⊗SM + SM⊗SP is the bond of `ring_model`.
"""
function cartesian_frame()
    ε(a, b, c) = (a - b) * (b - c) * (c - a) ÷ 2         # Levi-Civita on {1, 2, 3}
    S = [Complex{Int}[-im * ε(α, a, b) for a in 1:3, b in 1:3] for α in 1:3]
    W = Complex{Int}[-im 0 im; 1 0 1; 0 im 0]
    norms = [2, 1, 2]
    N = [i == j ? norms[i] : 0 for i in 1:3, j in 1:3]
    Splus, Sminus = S[1] + im * S[2], S[1] - im * S[2]
    @assert W' * W == N
    @assert W' * S[3] * W == N * SZ
    @assert W' * Splus * W == 2SP && W' * Sminus * W == 2SM
    @assert W' * S[3]^2 * W == N * [1 0 0; 0 0 0; 0 0 1]
    # S·S in the Cartesian frame equals Sᶻ⊗Sᶻ + (S⁺⊗S⁻ + S⁻⊗S⁺)/2 exactly.
    @assert 2 * sum(kron(s, s) for s in S) ==
            2kron(S[3], S[3]) + kron(Splus, Sminus) + kron(Sminus, Splus)
    return (; W, norms)
end

"Seed length (sites) of every stage-1 point."
const N0 = 6

"The stage-1 points: exact D and the seed inverse temperature β₀ (as in the Lean file)."
const POINTS = [(D = 4 // 1, beta0 = 9 // 4), (D = 1603 // 500, beta0 = 381 // 125)]

"""
Upper bounds on the seed defects (εₛ, εₜ) that `lean/CertifiedGaps/Models/SpinOneLargeD.lean`
assumes at each point (the first entries of its chain certificates `s` and `t`); the stage-1
certificates must not exceed them.
"""
const LEAN_SEED_BOUNDS = Dict(4 // 1 => (6411 // 10^5, 4053 // 10^5),
    1603 // 500 => (19164415 // 10^8, 46462453 // 10^9))

"`model_metadata` sources of the stage-1 certificates (written by `certify.jl`)."
const STAGE1_SOURCES = (dirs = (), files = ("src/model.jl", "certify.jl"))

"File stem of the stored results for an exact D, e.g. `D1603_500`."
stem(D) = "D$(numerator(D))_$(denominator(D))"

"Path (relative to this model directory) of the stored verified spectrum of the ring L at D."
spectrum_file(D, L) = "results/spectrum_$(stem(D))_L$L.json.gz"

"Path (relative to this model directory) of the stage-1 core certificate at D."
certificate_file(D) = "results/certificate_$(stem(D)).json"

end
