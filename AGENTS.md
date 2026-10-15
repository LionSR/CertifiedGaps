# Instructions for coding agents

This repository proves uniform spectral gaps of quantum spin chains with the spacetime
bootstrap (the two-purity doubling criterion; code identifiers keep the older "doubling" and
"two-purity" names) and certifies its finite inputs by computer. Read this file before
starting any task.

## The criterion (what a certificate must establish)

For a ring Hamiltonian H_L = Σ_j h_{j,j+1} (+ Σ_j o_j) with partition function Z_L(β):

    Pₛ(n, β) = Z_{2n}(β) / Z_n(β)²        (spatial purity)
    P_τ(n, β) = Z_n(2β) / Z_n(β)²          (thermal purity)

If (1) the bond is **reflection positive** in some local basis, i.e. −h = Σ_α C_α ⊗ conj(C_α)
with a real on-site term o (equivalently the realignment R[−h] is Hermitian PSD), and (2) the
**seed defects** εₛ ≥ 1 − Pₛ(n₀, β₀), εₜ ≥ 1 − P_τ(2n₀, β₀) at one even n₀ lie in the
admissible region (sufficient: both < 0.12731), then every even L ≥ n₀ has a unique ground
state and Δ_L ≥ `certified_gap_bound(εₛ, εₜ, β₀).bound`. Constant energy shifts cancel in both
purities. For c-site unit cells, lengths count cells and the cell may carry an internal
reflection F (−h = Σ C_α ⊗ F conj(C_α) F, o real and F-invariant; palindromic cell,
h_{c−k} = h_k), as in the scope note and the Lean `RPRing`. An odd n₀ needs an extra
argument (see `models/spin1_dimerized/`).

A **gap certificate** for a model consists of:

1. **RP proof**: explicit matrices C_α (exact rationals or algebraic numbers in closed form)
   with −h = Σ C_α ⊗ conj(C_α) verified by an exact identity, and o real. A floating-point
   `rp_test` is a screening step, not the proof.
2. **Seed enclosures**: rigorous upper bounds on 1 − Pₛ(n₀, β₀) and 1 − P_τ(2n₀, β₀), from
   rigorous enclosures of log Z_{n₀}(β₀), log Z_{2n₀}(β₀), log Z_{2n₀}(2β₀). Accepted sources:
   - exact rational/integer arithmetic, Arb balls (Nemo), or IntervalArithmetic;
   - verified eigenvalue enclosures of the full spectrum, or of every symmetry block with its
     multiplicity (summing to dⁿ): block entries enclosed exactly or in intervals, residual
     and orthogonality norms computed in `Float64` (BLAS allowed) **with an explicit a priori
     rounding-error bound** under the IEEE binary64 model (round-to-nearest, gradual
     underflow, the standard γₙ dot-product bound valid for any summation order and FMA),
     then a Weyl/polar-factor bound, and log Z summed in interval arithmetic. This is how
     every certificate outside `spin1/` was made (`CertifiedGaps.VerifiedSpectra` implements
     it); the README must state the arithmetic model as the trusted base.

   Not accepted: a `Float64` number without a proven error bound, an ARPACK/Lanczos run, DMRG,
   or an unverified dense eigensolve.
3. **The bound**: `CertifiedGaps.certified_gap_bound` (interval) or
   `CertifiedGaps.Recurrence.rational_gap_chain` (exact rational) applied to the seed bounds.
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
- **Smaller seeds.** The criterion needs only Pₛ(n₀, β₀) and P_τ(2n₀, β₀); by Theorem
  2.24(v) of the survey (kept locally, not in the repository), P_τ(2n₀) can be replaced by
  P_τ(n₀) and Pₛ(n₀, 2β₀), halving the largest ring. Unequal defects are cheap in εₛ and expensive in εₜ; spend effort accordingly.
- **Symmetry and exact counting.** Character formulas, SU(2) multiplet counting, sector
  inversion from twisted traces (as in `models/spin1/`) instead of enumerating states.
- **Variational sandwiches.** Upper bounds on E₀(2n) from trial states (MPS) and lower bounds
  on E₀(n) (e.g. Anderson/cluster bounds) bound the spatial defect δₙ = E₀(2n) − 2E₀(n).

Start every task by writing `models/<slug>/IDEAS.md`: the candidate routes, why each
might or might not work for this model, how its cost scales with n₀ and β₀, and which one you
chose. Commit that memo even if the implementation fails; a clear negative ("route X fails
because ...") is a result.

## Repository layout

- `outputs/`: LaTeX notes (main proof, scope, certificate method). Do not edit unless asked.
- `julia/`: shared core package `CertifiedGaps` (recurrence and certified bounds, exact
  arithmetic, results I/O, spin operators, RP tests and exact RP certificates, exact
  diagonalisation, verified ring spectra and `certify_gap`, Nemo/Arb certified MPO toolkit,
  ITensors purification extension). See `julia/README.md`.
- `models/spin1/`: exact-arithmetic reproduction of the spin-1 Heisenberg certificate (integer
  Poisson trace recursion with twisted rings and symmetry representatives, polynomial
  filters, MPS trial contractions). The template for large-n certificates.
- `models/<model-slug>/`: **where new certificates go** (one directory per model). The root
  `README.md` has the status table of all models; update its row when you add or change one.
- `lean/`: Lean 4 formalisation of the criterion and of each certified model's gap theorem
  (seeds are hypotheses).

## Directory format

Every `models/<slug>/` has exactly these top-level entries (examples:
`models/heisenberg_staggered_field` for a certificate, `models/blbq_haldane` for a negative
result):

| entry | content |
|---|---|
| `README.md` | required; sections below |
| `IDEAS.md` | required: the route memo written before computing; ends with `## Outcome` |
| `METHOD.md` | optional: the one long method/proof/audit document (no `REPORT.md` etc.) |
| `Project.toml`, `Manifest.toml` | the environment |
| `.gitignore` | local ignores for regenerable data (when there is any) |
| `certify.jl` | produces the certificate(s); omitted for a negative result |
| `check.jl` | verifies them: fast by default, `--full` for the slow replay |
| `src/` | library code (a module where natural) |
| `scripts/` | one-off, exploratory, scan, report and audit drivers |
| `test/runtests.jl` | fast tests (may include other files in `test/`) |
| `results/` | certificate data (`results/certificate*.json[.gz]`) and small summaries |

No other top-level files. The one exception is a file whose SHA-256 a certificate records by
path when that certificate cannot be regenerated in under 10 minutes: it stays where it is and
README "Layout" lists it with the reason. Never rename, move or edit a hashed file unless you
regenerate the certificate and every certified number is identical.

`README.md` uses these headings, in this order and with this wording:

1. `# <Model name in words>`: a 2–4 sentence summary, then a headline table
   (parameter | seed (n₀, β₀) | certified gap ≥ | lengths covered), or "No certificate" with
   the best numerical evidence.
2. `## Model`: Hamiltonian, conventions, cell structure, whether lengths count cells or sites.
3. `## Reflection positivity`: RP basis or frame, the exact identity and where it is checked.
4. `## Seeds and certified bounds`: exact (n₀, β₀, εₛ, εₜ) and the bound per parameter
   (`## Result` instead for a negative result).
5. `## Reproduce`: the certify, check, `check --full` and test commands with measured run times.
6. `## What is rigorous`: the trusted base (exact, interval, IEEE binary64 + BLAS rounding
   model), what is floating point only, platform caveats (e.g. LAPACK bit-identity).
7. `## Layout`: file map, including hash-pinned exceptions.
8. `## History` (optional): PRs, audits, superseded results.

Keep the README concise; long derivations go to `METHOD.md`.

## Environment

Julia 1.12. From the repository root:

```sh
julia --project=julia -e 'using Pkg; Pkg.instantiate(); Pkg.test()'
```

A new certificate directory is its own Julia project depending on the core:
`julia --project=models/<slug> -e 'using Pkg; Pkg.develop(path="julia")'` (adjust the
relative path to `../../julia` inside the project), or runs with `--project=julia`. Set
`OPENBLAS_NUM_THREADS=1` unless you need threaded BLAS.

## Rules

- Work only inside `models/<your-model-slug>/` unless the task says otherwise. If the
  core is missing something generic, add it to `julia/` in a separate, minimal commit with a
  test, and say why in the PR description.
- Certificates record SHA-256 hashes of `julia/src`, and some checkers compare them with the
  checkout. Any change to `julia/src` must keep every model checker passing; run them all.
- Never modify `models/spin1/results`, `models/spin2/results`, `models/aklt/results` or anything under `lean/`.
- No Python files: everything is Julia (and LaTeX for notes).
- Certificate JSON stores every exact number as a string (`"num//den"`); use
  `CertifiedGaps.Results.write_result` and `parse_exact_json`. Record the code's SHA-256
  (`source_hashes`) and the Julia/package versions (`run_metadata`).
- Every model directory follows "Directory format" above (files, `check.jl` modes, README
  headings).
- Report negative results too: if no admissible seed exists at accessible sizes, commit the
  scan with its (floating-point, clearly labelled) evidence and the estimated seed scale.
- Code style: 100 columns, docstrings on public functions, comments that say why.

## Open targets (floating point, from a prior survey; not certified)

Certified models are listed in the root `README.md`. Survey estimates for models without a
certificate:

| model | RP basis | seed (n₀, β₀) | float bound | exact gap |
|---|---|---|---|---|
| TFIM −ΣZZ − gΣX, g = 2 | σᶻ, C = Z | (128, 32) free fermions | 1.885 | 2 |
| TFIM g = 1.25 | σᶻ | (144, 64) | 0.451 | 0.5 |
| TFIM g = 1.1 | σᶻ | (384, 161) | 0.178 | 0.2 |
| cluster −ΣZXZ − gΣX, g = 0.5 (2-site cells) | cell basis CS(H⊗H) | (128, 64) cells | 0.932 | 1 |
| XY in a field, JₓJ_y < 0 | standard basis | free fermions | — | — |
| bilinear–biquadratic θ = −0.1 | Cartesian | ≈(34, 31) estimated | ≈0.107 | ≈0.338 |
