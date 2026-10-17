# The dimerized spin-1/2 Heisenberg ring as a core `RingModel`, and the model-specific stages
# built on the core's verified spectra: the all-β exclusion covers, the dimer perturbation
# bound and the parameter grid. Spectra, log Z enclosures, seed defects and the gap chain come
# from `CertifiedGaps` (`VerifiedSpectra`, `GapCertificates`).
module DimerizedHeisenberg

using LinearAlgebra, IntervalArithmetic
using CertifiedGaps.ExactArithmetic: RationalBig, round_down, parse_rational
using CertifiedGaps.RingModels: RingModel, rp_certificate, exact_matrix
using CertifiedGaps.VerifiedSpectra: BlockSpectrum, log_partition_enclosure,
    shifted_log_partition, ground_energy_enclosure
using CertifiedGaps.Recurrence: to_interval
using CertifiedGaps.Spins: exact_spin_half

export DELTAS, SEEDS, DIMER_BETA, rotated_model, rotated_rp, physical_terms, FRAME,
    frame_identity, MODEL_DIR, SOURCES, dimer_bounds, seed_defects, rejects, negative_cover,
    check_cover

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

const X, IY, Z = let S = exact_spin_half(Q); (S.x, S.iy, S.z) end   # IY = iSʸ, real
const ONE = Matrix{Q}(I, 2, 2)

"The swap of the two spins of a cell: the internal reflection F."
const SWAP = Q[1 0 0 0; 0 0 1 0; 0 1 0 0; 0 0 0 1]

"The spin flip X ⊗ X of a cell (up to a sign in the rotated frame): it reverses the charge."
const SPIN_FLIP = kron(2X, 2X)

"""
The frame rotation W = 1 ⊗ J of one cell, J = 2iSʸ = [0 1; −1 0] on its second spin: the
rotated ring Hamiltonian is U'H U with U = W^{⊗L}.
"""
const FRAME = kron(ONE, 2IY)

# S·S = SˣSˣ + SʸSʸ + SᶻSᶻ, with SʸSʸ = −(iSʸ)(iSʸ).
const HEIS = kron(X, X) - kron(IY, IY) + kron(Z, Z)

"""
    physical_terms(δ) -> (onsite, bond)

Physical cell term (1 + δ) S_a·S_b and bond (1 − δ) S_b·S_a' (second spin of a cell, first spin
of the next) on `kron(cell, cell)`.
"""
physical_terms(δ) = ((1 + δ) * HEIS, (1 - δ) * kron(ONE, HEIS, ONE))

# In the rotated frame JᵀSˣJ = −Sˣ, JᵀSᶻJ = −Sᶻ and Jᵀ(iSʸ)J = iSʸ.
rotated(J) = J * (-kron(X, X) - kron(IY, IY) - kron(Z, Z))

"""
    rotated_model(δ) -> RingModel

The ring of two-spin cells in the rotated frame: on-site term (1 + δ)(−SˣSˣ + SʸSʸ − SᶻSᶻ), bond
(1 − δ)(…) between the second spin of a cell and the first of the next. The cell charge
z_a − z_b is conserved and reversed by the spin flip of the cell.
"""
rotated_model(δ) = RingModel(4, kron(ONE, rotated(1 - δ), ONE), rotated(1 + δ);
    charge = [0, 1, -1, 0], flip = SPIN_FLIP,
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

Exact check that W = 1 ⊗ J is orthogonal and maps the physical cell term and bond onto those of
`rotated_model(δ)`. Since U = W^{⊗L} is a product of cell unitaries, U'H_L U is the rotated ring
Hamiltonian for every L, so both rings have the same spectrum.
"""
function frame_identity(δ)
    m = rotated_model(δ)
    onsite, bond = physical_terms(δ)
    W = FRAME
    WW = kron(W, W)
    return W' * W == I && exact_matrix(W' * onsite * W) == m.onsite &&
           exact_matrix(WW' * bond * WW) == m.bond
end

qinterval(x) = to_interval(x)

"""
    dimer_bounds(δ, n, β)

Exact norm-perturbation bounds about disconnected strong dimers (δ = 1): with
‖V_n‖ ≤ 3n(1 − δ)/4, every log Z error is ≤ β‖V_n‖. Enclosures of the defects they imply.
"""
function dimer_bounds(delta::Q, n::Int, beta::Q)
    v = 3n * (1 - delta) / 4
    x = exp(-qinterval(beta * (1 + delta)))
    single = (1 + 3x^2) / (1 + 3x)^2
    ps_lower = exp(-qinterval(4beta * v))
    pt_lower = single^(2n) * exp(-qinterval(8beta * v))
    return (; s = 1 - ps_lower, t = 1 - pt_lower, norm_small = v)
end

"""
    seed_defects(small, large, β) -> (; s, t)

Enclosures of 1 − Pₛ(n, β) and 1 − P_τ(2n, β) from the verified spectra of the rings n and 2n
(core `log_partition_enclosure`) at an exact β.
"""
function seed_defects(small, large, β)
    zn, z2n = log_partition_enclosure(small, β), log_partition_enclosure(large, β)
    z2n2 = log_partition_enclosure(large, 2β)
    return (; s = 1 - exp(z2n - 2zn), t = 1 - exp(z2n2 - 2z2n))
end

# Any exact shift is valid; the smallest center keeps the exponents near or below zero.
shift(spectrum) = minimum(minimum(s.centers) for s in spectrum)

# Rational lower bounds on both defects on the closed box [a, b] of β. The energy shifts of the
# two rings are combined exactly, so the width of the β box does not multiply E₀.
function lower_defects(small, large, a::Q, b::Q)
    β = interval(inf(qinterval(a)), sup(qinterval(b)))
    sn, s2n = shift(small), shift(large)
    zn = shifted_log_partition(small, β, sn)
    z2n = shifted_log_partition(large, β, s2n)
    z2n2 = shifted_log_partition(large, 2β, s2n)
    ps = exp(-β * qinterval(Q(s2n) - 2Q(sn)) + z2n - 2zn)
    pt = exp(z2n2 - 2z2n)
    return max(zero(Q), Q(inf(1 - ps))), max(zero(Q), Q(inf(1 - pt)))
end

"""
    rejects(s, t; steps = 16)

Doubling step at which lower bounds (s, t) on the defects provably leave the admissible region
(z ≥ 1/2), or `nothing`. The map is monotone, so rounding down after every step keeps lower
bounds.
"""
function rejects(s::Q, t::Q; steps = 16)
    grid = big(10)^15
    for j in 0:steps
        z = 1 - (1 - s)^2 * (1 - t)
        z >= 1 // 2 && return j
        s = round_down(z^2 / (2 * (1 - z)^2), grid)
        zt = 1 - (1 - t)^2 * (1 - s)
        zt >= 1 // 2 && return j + 1  # The next z(s', t') ≥ zt²/(2(1 − zt)²) ≥ 1/2.
        t = round_down(zt^2 / (2 * (1 - zt)^2), grid)
    end
    return nothing
end

# Rigorous lower bound on d = E₀(2n) − 2E₀(n) from the core ground-energy enclosures.
energy_defect_lower(small, large) =
    ground_energy_enclosure(large)[1] - 2ground_energy_enclosure(small)[2]

# Smallest dyadic B ≥ 32 with 2^{4n} e^{−B d} ≤ 1/2, so Pₛ(n, β) ≤ 1/2 for all β ≥ B.
function tail_holds(n, B::Q, d::Q)
    d > 0 && sup(exp(qinterval(4n) * log(interval(2)) - qinterval(B * d))) <= 1 // 2
end

"""
    negative_cover(small, large, n; max_depth = 16) -> Dict

Proof that every β₀ > 0 fails the doubling recurrence at the even seed length n, from the
verified spectra of the rings n and 2n:
(0, 1] by monotonicity of P_τ(2n, β) once P_τ(2n, 1) ≤ 1/2; [B, ∞) by
Pₛ(n, β) ≤ 2^{4n} e^{−β d}, d = E₀(2n) − 2E₀(n) > 0; and [1, B] by boxes on which lower bounds
on both defects are rejected by `rejects`.
"""
function negative_cover(small, large, n; max_depth = 16)
    first = seed_defects(small, large, 1 // big(1))
    inf(first.t) >= 1 // 2 || return Dict("status" => "unresolved_low_beta")
    d = energy_defect_lower(small, large)
    d > 0 || return Dict("status" => "unresolved_high_beta")
    tail = big(32) // 1
    while !tail_holds(n, tail, d)
        tail *= 2
    end
    cover, unresolved = Any[], Any[]
    function visit(a, b, depth)
        s, t = lower_defects(small, large, a, b)
        j = rejects(s, t)
        if j !== nothing
            push!(cover, Dict("a" => string(a), "b" => string(b), "s_lower" => string(s),
                "t_lower" => string(t), "rejected_by_step" => string(j)))
        elseif depth == max_depth
            push!(unresolved, [string(a), string(b)])
        else
            mid = (a + b) / 2
            visit(a, mid, depth + 1)
            visit(mid, b, depth + 1)
        end
    end
    visit(big(1) // 1, tail, 0)
    return Dict("status" => isempty(unresolved) ? "all_beta_rejected" : "unresolved",
        "n0_cells" => string(n), "low_beta_end" => "1//1",
        "thermal_defect_lower_at_beta_1" => string(Q(inf(first.t))), "tail_beta" => string(tail),
        "ground_energy_defect_lower" => string(d), "cover" => cover, "unresolved" => unresolved,
        "scope" => "No positive rate from the stated doubling recurrence, not a gaplessness claim")
end

"""
    check_cover(record; small = nothing, large = nothing)

Check an all-β exclusion. Without spectra: the stored lower bounds (P_τ(2n, 1) ≤ 1/2, the
energy defect and its tail, every box's defects) reject the recurrence and the boxes tile
[1, B]. With the verified spectra of the rings n and 2n: every bound is recomputed from them
and the stored values are not used.
"""
function check_cover(record; small = nothing, large = nothing)
    fresh = small !== nothing
    record["status"] == "all_beta_rejected" || error("unresolved beta region")
    isempty(record["unresolved"]) || error("unresolved boxes")
    n = parse(Int, string(record["n0_cells"]))
    t1 = fresh ? Q(inf(seed_defects(small, large, 1 // big(1)).t)) :
         parse_q(record["thermal_defect_lower_at_beta_1"])
    t1 >= 1 // 2 || error("small-beta exclusion failed")
    d = fresh ? energy_defect_lower(small, large) :
        parse_q(record["ground_energy_defect_lower"])
    tail = parse_q(record["tail_beta"])
    tail_holds(n, tail, d) || error("large-beta exclusion failed")
    previous = 1 // big(1)
    for box in record["cover"]
        a, b = parse_q(box["a"]), parse_q(box["b"])
        a == previous < b || error("beta cover has a gap or overlap")
        s, t = fresh ? lower_defects(small, large, a, b) :
               (parse_q(box["s_lower"]), parse_q(box["t_lower"]))
        rejects(s, t) !== nothing || error("beta box is not excluded")
        previous = b
    end
    previous == tail || error("beta cover misses its tail")
    return true
end

parse_q(x) = parse_rational(x)

const MODEL_DIR = normpath(joinpath(@__DIR__, ".."))

"""
The model-local sources hashed into every result by `Results.model_metadata` and enforced by
`Results.verify_model_metadata`: `src/`, `certify.jl`, and `check.jl` (the environment is the
core's, recorded under `core_sha256`).
"""
const SOURCES = (; dirs = ("src",), files = ("certify.jl", "check.jl"))

end
