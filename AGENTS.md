# Instructions for coding agents

This repository proves uniform spectral gaps of quantum spin chains using the **spacetime
bootstrap** (also called the *two-purity doubling criterion*; code identifiers retain the shorter
`doubling` and `two_purity` names) and certifies its finite numerical inputs by computer. Read this
file before starting any task.

## The criterion (what a certificate must establish)

For a periodic ring Hamiltonian $H_L = \sum_{j=1}^L h_{j,j+1} + \sum_{j=1}^L o_j$ of $L$ unit cells
with partition function $Z_L(\beta) = \operatorname{Tr} e^{-\beta H_L}$, define the spatial and thermal
purities:

    Pₛ(n, β) = Z_{2n}(β) / Z_n(β)²        (spatial purity)
    P_τ(n, β) = Z_n(2β) / Z_n(β)²          (thermal purity)

If:
1. the nearest-neighbor bond interaction $h$ is **reflection positive** in some local basis, meaning
   $-h = \sum_\alpha C_\alpha \otimes \overline{C_\alpha}$ with a real symmetric on-site term $o$
   (equivalently, the bond realignment $R[-h]$ is Hermitian positive semidefinite), and
2. the **seed purity defects** $\varepsilon_{\mathrm{s}} \ge 1 - P_{\mathrm{s}}(n_0, \beta_0)$ and
   $\varepsilon_\tau \ge 1 - P_\tau(2n_0, \beta_0)$ at one even seed length $n_0$ lie in the admissible
   region of the doubling recurrence (a sufficient condition is $\max(\varepsilon_{\mathrm{s}}, \varepsilon_\tau) < \varepsilon_* \approx 0.12731$),

then every even ring length $L \ge n_0$ has a unique ground state and a spectral gap
$\Delta_L = E_1(L) - E_0(L) \ge$ `certified_gap_bound(εₛ, εₜ, β₀).bound`. Constant energy shifts cancel
identically in both purities. For $c$-site unit cells, ring lengths $L$ count unit cells and the cell
may carry an internal reflection $F$ (a real symmetric involution satisfying $-h = \sum_\alpha C_\alpha \otimes F \overline{C_\alpha} F$ and $F o F = o$, corresponding to a palindromic cell with $h_{c-k} = h_k$), as proved in `notes/thermal-mpo-gap-scope.tex` and formalized in the Lean structure `RPRing`. Using an odd seed length $n_0$ requires an explicit positivity argument for the transfer operator (see `models/spin1_dimerized/`).

A **gap certificate** for a model consists of four components:

1. **Reflection-positivity proof**: explicit matrices $C_\alpha$ (with exact rational or closed-form
   algebraic entries) satisfying $-h = \sum_\alpha C_\alpha \otimes \overline{C_\alpha}$ (or $-h = \sum_\alpha C_\alpha \otimes F \overline{C_\alpha} F$ with an internal cell reflection $F$) verified by an exact algebraic identity, together with a real symmetric on-site term $o$. A floating-point `rp_test` is a screening step, not a proof.
2. **Seed enclosures**: rigorous rational upper bounds $\varepsilon_{\mathrm{s}} \ge 1 - P_{\mathrm{s}}(n_0, \beta_0)$
   and $\varepsilon_\tau \ge 1 - P_\tau(2n_0, \beta_0)$, obtained from rigorous enclosures of
   $\log Z_{n_0}(\beta_0)$, $\log Z_{2n_0}(\beta_0)$, and $\log Z_{2n_0}(2\beta_0)$. Accepted sources are:
   - exact rational or integer arithmetic, Arb ball arithmetic (Nemo), or `IntervalArithmetic`;
   - verified eigenvalue enclosures of the full spectrum, or of every symmetry block with its exact
     multiplicity summing to $d^n$: block entries enclosed exactly or in intervals, residual and
     orthogonality norms computed in `Float64` (BLAS allowed) **with an explicit a priori
     rounding-error bound** under the IEEE binary64 model (round-to-nearest, gradual underflow, and
     the standard $\gamma_n$ dot-product bound valid for any summation order and FMA), followed by a
     Weyl/Ostrowski/polar-factor perturbation bound and $\log Z$ summation in interval arithmetic.
     `CertifiedGaps.VerifiedSpectra` implements this pipeline; the model's `README.md` must state this
     arithmetic model as its trusted base.

   *Not accepted as a certificate:* a `Float64` number without a proven error bound, an ARPACK/Lanczos
   run, DMRG, thermal MPO purification without rigorous truncation bounds, or an unverified dense eigensolve.
3. **The gap bound**: `CertifiedGaps.certified_gap_bound` (interval arithmetic) or
   `CertifiedGaps.Recurrence.rational_gap_chain` (exact rational certificate) applied to $(\varepsilon_{\mathrm{s}}, \varepsilon_\tau, \beta_0)$.
4. **A standalone checker** (`check.jl`) that recomputes steps 1–3 from the stored certificate in a few
   minutes (and supports `--full` for a full spectral replay) and fails loudly on any mismatch.

## Terminology and notation conventions

Use consistent terminology and notation across all documentation and comments:

- **Method name:** *spacetime bootstrap* or *two-purity doubling criterion*.
- **Purities and defects:** *spatial purity* $P_{\mathrm{s}}(n, \beta) = Z_{2n}(\beta)/Z_n(\beta)^2$ (`Pₛ(n, β)`)
  and *thermal purity* $P_\tau(n, \beta) = Z_n(2\beta)/Z_n(\beta)^2$ (`P_τ(n, β)`; never "Gibbs purity" or
  "temporal purity"). Their complements $1 - P_{\mathrm{s}}$ and $1 - P_\tau$ are the *spatial defect* and
  *thermal defect*, bounded at the seed scale $(n_0, \beta_0)$ by $(\varepsilon_{\mathrm{s}}, \varepsilon_\tau)$ (`εₛ`, `εₜ`).
- **Hamiltonian and transfer operator:** $H_L = \sum_j h_{j,j+1} + \sum_j o_j$ on a ring of $L$ unit cells
  ($N = cL$ spins for a $c$-site cell), with bond term $h$, on-site term $o$, internal cell reflection $F$,
  and spatial transfer operator $X_\beta$ ($\operatorname{Tr} X_\beta^L = Z_L(\beta)$). For spin-$S$ chains in the
  Cartesian basis, write $G_\alpha = i S^\alpha$ ($\alpha \in \{x, y, z\}$) for the real antisymmetric spin
  generators satisfying $-\mathbf{S}_j \cdot \mathbf{S}_{j+1} = \sum_\alpha G_{\alpha,j} \otimes G_{\alpha,j+1}$.
- **Spelling:** Use standard American spelling (`diagonalization`, `dimerized`, `normalization`,
  `magnetization`, `formalized`, `optimized`, `parameter`, `neighbor`, `labeled`) throughout Markdown prose.

## Approach: think before computing

Brute-force exact diagonalization stops near 20 spin-½ sites, whereas sharp seeds require rings of
roughly $20\xi$ sites ($\xi$ the correlation length) at $\beta_0 \approx (5\text{–}11)/\Delta$. A task that only
diagonalizes larger rings is a baseline, not a scalable method. Look first for mathematical structure
that makes the seed purities computable or boundable beyond raw dense diagonalization:

- **Closed forms:** Free fermions, commuting or frustration-free points, and products of small blocks
  give $Z_L(\beta)$ in closed form at arbitrary $L$.
- **Expansions with rigorous remainders:** Convergent cluster or polymer expansions (Kotecký–Preiss),
  high- or low-temperature expansions, and perturbations around solvable limits yield
  $\log Z_L(\beta) = L f(\beta) + O(e^{-L/\xi})$ on rings, so $P_{\mathrm{s}} = \exp(\log Z_{2n} - 2\log Z_n)$ is $1$
  up to an exponentially small, rigorously bounded remainder uniformly in $n$.
- **Transfer-operator structure:** For a self-adjoint transfer operator $X_\beta$, the spatial defect
  is bounded by $\sum_{i \ge 1} (\lambda_i/\lambda_0)^n$. Any rigorous bound on the transfer-operator
  spectral ratio (via Perron–Frobenius and Birkhoff contraction for positive kernels, kernel comparison,
  or rank-one perturbation theory) bounds $P_{\mathrm{s}}$ for all $n$ simultaneously.
- **Perturbation from a solvable point:** Start where the seed is exact (dimer limit, large field,
  strong rungs, large single-ion anisotropy $D$) and propagate a certified error bound (Duhamel,
  Peierls–Bogoliubov, Weyl, or operator inequalities) to the target coupling while tracking the
  admissible defect region ($< 0.12731$).
- **Asymmetric and half-length seeds:** The criterion needs only $P_{\mathrm{s}}(n_0, \beta_0)$ and
  $P_\tau(2n_0, \beta_0)$. By Theorem 2.24(v) of the local survey notes, $P_\tau(2n_0, \beta_0)$ can be
  bounded in terms of $P_\tau(n_0, \beta_0)$ and $P_{\mathrm{s}}(n_0, 2\beta_0)$, halving the largest required
  ring length. Moreover, the asymmetric recurrence is much more tolerant of the spatial defect
  $\varepsilon_{\mathrm{s}}$ than the thermal defect $\varepsilon_\tau$; allocate computational effort accordingly.
- **Symmetry and exact counting:** Use character formulas, $\mathrm{SU}(2)$ multiplet counting, and
  sector inversion from twisted traces (as in `models/spin1/`) rather than enumerating individual basis states.
- **Variational sandwiches:** Upper bounds on $E_0(2n)$ from trial states (such as MPS) combined with
  lower bounds on $E_0(n)$ (such as Anderson cluster bounds) bound the ground-energy doubling mismatch
  $\delta(n) = E_0(2n) - 2E_0(n)$.

Start every task by writing `models/<slug>/IDEAS.md`: explain the candidate routes, why each might or
might not work for the model, how its cost scales with $n_0$ and $\beta_0$, and which route was chosen.
Commit that memo even if the implementation fails; a clear negative result ("route X fails because ...")
is valuable documentation.

## Repository layout

- `notes/`: LaTeX notes (main proof, scope, certificate method, and $2+1$D extension). Do not edit unless asked.
- `julia/`: shared core package `CertifiedGaps` (recurrence and certified bounds, exact arithmetic,
  results I/O, spin operators, RP tests and exact RP certificates, exact diagonalization, verified ring
  spectra and `certify_gap`, Nemo/Arb certified MPO toolkit, and ITensors purification extension). See
  `julia/README.md`.
- `models/spin1/`: exact-arithmetic reproduction of the spin-1 Heisenberg certificate (integer Poisson
  trace recursion with twisted rings and symmetry representatives, polynomial filters, and MPS trial
  contractions). Serves as the reference template for large-$n$ exact certificates.
- `models/<model-slug>/`: **where new model certificates live** (one directory per model). Update the
  summary table in the root `README.md` whenever you add or update a model.
- `lean/`: Lean 4 formalization of the criterion and of each certified model's gap theorem (library
  `CertifiedGaps`, with the two seeds as hypotheses), plus the opt-in library `SeedKernel` for
  kernel-checked seed arithmetic (`lake build SeedKernel`).
- `models/lean_seeds/`: generators and studies for seeds checked inside the Lean kernel (`decide +kernel`
  only; Julia generates witness files and Lean verifies them); its generators write the witness files
  under `lean/SeedKernel`.

## Directory format

Every `models/<slug>/` directory must contain only these top-level entries (see
`models/heisenberg_staggered_field` for a standard certificate and `models/blbq_haldane` for a negative
result):

| entry | content |
|---|---|
| `README.md` | required; must use the section headings listed below |
| `IDEAS.md` | required: the route analysis memo written before computing, ending with `## Outcome` |
| `METHOD.md` | optional: the single detailed method, proof, and audit document (do not create separate `REPORT.md`, `NOTES.md`, etc.) |
| `Project.toml`, `Manifest.toml` | only for a model that needs packages outside the core environment, or whose immutable results pin its own environment (`spin1`, `spin2`, `aklt`, `lean_seeds`); every other model runs in the core environment (see "Environment") |
| `.gitignore` | local ignores for regenerable data (if any) |
| `certify.jl` | generates the certificate(s); omitted for a negative result |
| `check.jl` | verifies the certificate(s): fast rational/hash verification by default, `--full` for the full spectral replay |
| `src/` | model-specific library code |
| `scripts/` | exploratory, scan, report, and audit driver scripts |
| `test/runtests.jl` | fast unit and regression tests (may include auxiliary files in `test/`) |
| `results/` | committed certificate data (`results/certificate*.json[.gz]`) and compact summaries |

Do not add other top-level files. The only exception is a file whose SHA-256 hash is recorded by path
inside a committed certificate that takes more than 10 minutes to regenerate: keep such files in place
and document them in the `## Layout` section of `README.md`. Never rename, move, or edit a hash-pinned
file unless you regenerate the certificate and verify that every certified number is identical.

Every model `README.md` must use these headings, in this order and with this exact wording:

1. `# <Model name in words>`: a 2–4 sentence summary, followed by a headline table
   (`parameter | seed (n₀, β₀) | certified gap ≥ | lengths covered`), or `"No certificate"` with a
   summary of the best numerical evidence.
2. `## Model`: the Hamiltonian, conventions, unit-cell structure, and whether lengths $L$ count cells or sites.
3. `## Reflection positivity`: the reflection-positive basis or rotated frame, the exact factorization
   identity, and where it is verified.
4. `## Seeds and certified bounds`: the exact values of $(n_0, \beta_0, \varepsilon_{\mathrm{s}}, \varepsilon_\tau)$ and the
   resulting gap bound per parameter point (use `## Result` instead for a negative result or pilot study).
5. `## Reproduce`: the `certify.jl`, `check.jl`, `check.jl --full`, and `test/runtests.jl` commands with
   measured run times.
6. `## What is rigorous`: the trusted base (exact arithmetic, interval arithmetic, IEEE binary64 + BLAS
   rounding model), which scripts or diagnostics are floating-point only, and any platform caveats.
7. `## Layout`: file map of the directory, noting any hash-pinned files.
8. `## History` (optional): summary of PRs, audits, or superseded runs.

Keep `README.md` concise and self-contained; place long derivations, scan tables, and audit logs in `METHOD.md`.

## Environment

Use Julia 1.12. From the repository root:

```sh
julia --project=julia -e 'using Pkg; Pkg.instantiate(); Pkg.test()'
```

Models share one environment, the core's: run every model script from the repository root with
`--project=julia` (for example `julia --project=julia models/<slug>/check.jl`). A model does not
get its own `Project.toml`/`Manifest.toml`; if it needs a package the core lacks, add the
dependency to the core in a separate commit, or explain in its `README.md` why it cannot share the
core environment. Set `OPENBLAS_NUM_THREADS=1` unless multi-threaded BLAS is specifically needed.

A gap certificate for a model that is a `RingModel` goes through the shared pipeline:
`GapCertificates.certify_gap` (exact RP check, `VerifiedSpectra` block spectra with
`verified_eigen` enclosures, seed defects, exact chain) and `check_certificate` in `check.jl`.
Symmetry reductions (charge, `flip`, on-site `symmetries`, cell `blocking`), reuse and storage of
spectra (`seed_bounds(...; spectra)`, `write_spectrum`), and gap enclosures
(`spectral_gap_enclosure`) are core features; do not reimplement them in a model.

## Rules

- Work only inside `models/<your-model-slug>/` unless the task explicitly specifies otherwise. If
  `CertifiedGaps` is missing generic functionality needed by multiple models, add it to `julia/` in a
  separate, minimal commit with unit tests and explain why.
- **Source hashes:** a checker enforces (fails on mismatch) only **model-local** sources, that is,
  the model's own `src/`, `certify.jl`, `check.jl`, the scripts that produce the certificate, and
  its `Project.toml`/`Manifest.toml` if it has its own environment, recorded in
  `metadata["source_sha256"]`; if a checker checks
  the recorded key set, it compares it with that model-local set only. The shared core
  (`julia/src`, `julia/Project.toml`, `julia/Manifest.toml`) is recorded under
  `metadata["core_sha256"]`, keyed by repository path (`julia/src/Results.jl`), as provenance
  only: the checker reruns the core arithmetic instead of pinning its bytes. Do not hash
  `notes/` or `METHOD.md`. Editing or adding core files, notes, or a `METHOD.md` therefore does
  not invalidate a certificate, but checkers rerun the core arithmetic, so any change to
  `julia/src` must keep every model checker and test passing; run them all. Use
  `Results.model_metadata` to record both hash sets and `Results.verify_model_metadata` (same
  `dirs` and `files`) to check the model-local set.
- Never modify `models/spin1/results`, `models/spin2/results`, `models/aklt/results`, or anything under
  `lean/` (except when running `models/lean_seeds/certify.jl` to regenerate the witness files in
  `lean/SeedKernel`).
- Do not add Python files; all code in this repository is written in Julia, Lean 4, or LaTeX.
- Store every exact rational or integer in certificate JSON files as a string (`"num//den"` or `"int"`)
  using `CertifiedGaps.Results.write_result` and `parse_exact_json`. Record source file SHA-256 hashes
  (`source_hashes`, split as in "Source hashes" above) and Julia/package versions (`run_metadata`).
- Every model directory follows "Directory format" above (files, `check.jl` modes, and `README.md`
  headings).
- Report negative results clearly: if no admissible seed exists at accessible system sizes, commit the
  scan scripts and summary data (clearly labeled as floating-point diagnostics) together with the
  estimated seed scale $(n_0, \beta_0)$ in `README.md` and `IDEAS.md`.
- Code style: 100-column line limit, docstrings on public functions, and comments that explain *why*
  rather than restating *what* the code does.

## Open targets (floating-point survey estimates; not certified)

Certified models are listed in the root `README.md`. The table below gives floating-point seed estimates
from a prior survey for models that do not yet have a certificate:

| model | RP basis | estimated seed $(n_0, \beta_0)$ | float bound | exact gap |
|---|---|---|---|---|
| TFIM $-\sum ZZ - g\sum X$, $g = 2$ | $\sigma^z$, $C = Z$ | $(128, 32)$ free fermions | $1.885$ | $2$ |
| TFIM $g = 1.25$ | $\sigma^z$ | $(144, 64)$ | $0.451$ | $0.5$ |
| TFIM $g = 1.1$ | $\sigma^z$ | $(384, 161)$ | $0.178$ | $0.2$ |
| cluster $-\sum ZXZ - g\sum X$, $g = 0.5$ (2-site cells) | cell basis $\mathrm{CS}(H \otimes H)$ | $(128, 64)$ cells | $0.932$ | $1$ |
| XY chain in a field, $J_x J_y < 0$ | standard basis | free fermions | — | — |
| spin-1 bilinear–biquadratic $\theta = -0.1$ | Cartesian | $\approx (34, 31)$ estimated | $\approx 0.107$ | $\approx 0.338$ |
