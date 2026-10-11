# Instructions for coding agents

This repository proves uniform spectral gaps of quantum spin chains with the two-purity
doubling criterion, and certifies the finite inputs with exact or outward-rounded arithmetic.
Read this file before starting any task.

## The criterion (what a certificate must establish)

For a ring Hamiltonian H_L = Σ_j h_{j,j+1} (+ Σ_j o_j) with partition function Z_L(β):

    Pₛ(n, β) = Z_{2n}(β) / Z_n(β)²        (spatial purity)
    P_τ(n, β) = Z_n(2β) / Z_n(β)²          (thermal purity)

If (1) the bond is **reflection positive** in some local basis, i.e. −h = Σ_α C_α ⊗ conj(C_α)
with a real on-site term o (equivalently the realignment R[−h] is Hermitian PSD), and (2) the
**seed defects** εₛ ≥ 1 − Pₛ(n₀, β₀), εₜ ≥ 1 − P_τ(2n₀, β₀) at one even n₀ lie in the
admissible region (sufficient: both < 0.12731), then every even L ≥ n₀ has a unique ground
state and Δ_L ≥ `certified_gap_bound(εₛ, εₜ, β₀).bound`. Constant energy shifts cancel in both
purities. For c-site unit cells, lengths count cells and the cell must satisfy the reflection
condition of `notes`/scope note (palindromic cell, h_{c−k} = h_k).

A **gap certificate** for a model consists of:

1. **RP proof**: explicit matrices C_α (exact rationals or algebraic numbers in closed form)
   with −h = Σ C_α ⊗ conj(C_α) verified by an exact identity, and o real. A floating-point
   `rp_test` is a screening step, not the proof.
2. **Seed enclosures**: rigorous upper bounds on 1 − Pₛ(n₀, β₀) and 1 − P_τ(2n₀, β₀), from
   rigorous enclosures of log Z_{n₀}(β₀), log Z_{2n₀}(β₀), log Z_{2n₀}(2β₀). Acceptable sources:
   exact rational/integer arithmetic, Arb balls (Nemo), IntervalArithmetic, or verified
   eigenvalue enclosures (e.g. Weyl/residual bounds computed in interval arithmetic) of the
   full spectrum or of every symmetry block with its multiplicity. **No inequality may be
   accepted on the strength of a `Float64` computation**, an ARPACK/Lanczos run, DMRG, or an
   unverified dense eigensolve.
3. **The bound**: `ThermalMPOGap.certified_gap_bound` (interval) or
   `ThermalMPOGap.Recurrence.rational_gap_chain` (exact rational) applied to the seed bounds.
4. **A checker** that recomputes 1–3 from the stored certificate in a few minutes and fails
   loudly on any mismatch.

## Approach: think before computing

Brute force does not scale. Exact diagonalisation stops near 20 spins, while sharp seeds
need rings of about 20ξ sites (ξ the correlation length) at β₀ ≈ (5–11)/Δ. A task that only
diagonalises larger rings is a baseline, not a result. Look first for structure that makes the
seed purities computable or boundable at sizes far beyond ED:

- **Closed forms.** Free fermions, commuting or frustration-free points, products of small
  blocks: Z_L(β) exactly, at any L.
- **Expansions with rigorous remainders.** Convergent cluster or polymer expansions
  (Kotecký–Preiss), high- or low-temperature expansions, Kennedy–Tasaki-type expansions around
  solvable points. These give log Z_L(β) = L f(β) + O(e^{−L/ξ}) on rings, so the spatial
  purity Pₛ = exp(log Z₂ₙ − 2 log Zₙ) is 1 up to an exponentially small, *bounded* error,
  uniformly in n.
- **Transfer-operator structure.** The spatial defect is Σ_{i≥1} (λᵢ/λ₀)ⁿ for the
  self-adjoint transfer operator X_β. Any rigorous bound on its gap (Perron–Frobenius and
  Birkhoff contraction for positive kernels, comparison with a solvable kernel, perturbation
  of a rank-one kernel) bounds Pₛ for all n at once.
- **Perturbation from a solvable point.** Start where the seed is exact (dimer limit, large
  field, strong rungs, large D) and carry a certified error (Duhamel, Peierls–Bogoliubov,
  Weyl, operator inequalities) to the model of interest. Track how far the admissible region
  (defects < 0.1273) lets you go.
- **Smaller seeds.** The criterion needs only Pₛ(n₀, β₀) and P_τ(2n₀, β₀); by survey
  Theorem 2.24(v), P_τ(2n₀) can be replaced by P_τ(n₀) and Pₛ(n₀, 2β₀), halving the largest
  ring. Unequal defects are cheap in εₛ and expensive in εₜ; spend effort accordingly.
- **Symmetry and exact counting.** Character formulas, SU(2) multiplet counting, sector
  inversion from twisted traces (as in `spin1/`) instead of enumerating states.
- **Variational sandwiches.** Upper bounds on E₀(2n) from trial states (MPS) and lower bounds
  on E₀(n) (e.g. Anderson/cluster bounds) bound the spatial defect δₙ = E₀(2n) − 2E₀(n).

Start every task by writing `certificates/<slug>/IDEAS.md`: the candidate routes, why each
might or might not work for this model, how its cost scales with n₀ and β₀, and which one you
chose. Commit that memo even if the implementation fails; a clear negative ("route X fails
because ...") is a result.

## Repository layout

- `outputs/`: LaTeX notes (main proof, scope, certificate method). Do not edit unless asked.
- `julia/`: shared core package `ThermalMPOGap` (recurrence and certified bounds, exact
  arithmetic, results I/O, spin operators, RP test, exact diagonalisation, models, Nemo/Arb
  certified MPO toolkit, ITensors purification extension). See `julia/README.md`.
- `spin1/`: exact-arithmetic reproduction of the spin-1 Heisenberg certificate (integer
  Poisson trace recursion with twisted rings and symmetry representatives, polynomial
  filters, MPS trial contractions). The template for large-n certificates.
- `spin2/`, `aklt/`: spin-2 and 2D AKLT investigations (no certificate yet).
- `lean/`: Lean 4 formalisation of the abstract criterion (seeds are hypotheses).
- `certificates/<model-slug>/`: **where new certificates go** (one directory per model).

## Environment

Julia 1.12. From the repository root:

```sh
julia --project=julia -e 'using Pkg; Pkg.instantiate(); Pkg.test()'
```

A new certificate directory is its own Julia project depending on the core:
`julia --project=certificates/<slug> -e 'using Pkg; Pkg.develop(path="julia")'` (adjust the
relative path to `../../julia` inside the project). Set `OPENBLAS_NUM_THREADS=1` unless you
need threaded BLAS.

## Rules

- Work only inside `certificates/<your-model-slug>/` unless the task says otherwise. If the
  core is missing something generic, add it to `julia/` in a separate, minimal commit with a
  test, and say why in the PR description.
- Never modify `spin1/results`, `spin2/results`, `aklt/results` or anything under `lean/`.
- No Python files: everything is Julia (and LaTeX for notes).
- Certificate JSON stores every exact number as a string (`"num//den"`); use
  `ThermalMPOGap.Results.write_result` and `parse_exact_json`. Record the code's SHA-256
  (`source_hashes`) and the Julia/package versions (`run_metadata`).
- Each certificate directory has: `Project.toml`, `src/` or a single `certify.jl` that produces
  `certificate.json`, `check.jl` that verifies it, `test/runtests.jl` (fast), and a `README.md`
  stating the model, the RP basis, the seed (n₀, β₀), the certified bound, run time, and what
  is and is not rigorous.
- Report negative results too: if no admissible seed exists at accessible sizes, commit the
  scan with its (floating-point, clearly labelled) evidence and the estimated seed scale.
- Code style: 100 columns, docstrings on public functions, comments that say why.

## Known numbers (floating point, from a prior survey; to be certified)

| model | RP basis | seed (n₀, β₀) | float bound | exact gap |
|---|---|---|---|---|
| TFIM −ΣZZ − gΣX, g = 2 | σᶻ, C = Z | (128, 32) free fermions | 1.885 | 2 |
| TFIM g = 1.25 | σᶻ | (144, 64) | 0.451 | 0.5 |
| TFIM g = 1.1 | σᶻ | (384, 161) | 0.178 | 0.2 |
| cluster −ΣZXZ − gΣX, g = 0.5 (2-site cells) | cell basis CS(H⊗H) | (128, 64) cells | 0.932 | 1 |
| Heisenberg + staggered field h_s = 1 | sublattice-rotated, C = {Sˣ, iSʸ, Sᶻ} | (8, 5.25) | 1.149 | ≈1.874 |
| dimerised Heisenberg δ = 0.5 (2-site cells) | rotated frame | (4, 6.25) cells | 0.533 | ≈1.19 |
| dimerised Heisenberg δ = 0.3 | rotated frame | (4, 6.5) cells | 0.172 | ≈0.81 |
| spin-1 Heisenberg + D(Sᶻ)², D = 4 | Cartesian, o = D·diag(1,1,0) | (6, 3) | 0.544 | ≈2.27 |
| two-leg ladder J⊥ = 4 (1-rung cells) | checkerboard sublattice rotation | (4, 2) rungs | 1.102 | ≈3.14 |
| XY in a field, JₓJ_y < 0 | standard basis | free fermions | — | — |
| bilinear–biquadratic θ = −0.1 | Cartesian | ≈(34, 31) estimated | ≈0.107 | ≈0.338 |
| spin-1 Heisenberg (the note) | Cartesian | (2304, 784) certified in `spin1/`; gap > 0.0729 for L ≥ 36 claimed in `spin1/experiments` | | 0.4105 |
