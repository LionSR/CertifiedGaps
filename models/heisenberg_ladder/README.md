# Two-leg spin-1/2 Heisenberg ladder

This module certifies uniform spectral gaps for the two-leg spin-1/2 Heisenberg ladder at rung couplings $J_\perp = 4$, $7/2$, $3$, $5/2$, and $163069/65536$ from verified complete spectra of rings of $2$, $4$, and $8$ rungs. Each certificate establishes a unique ground state and $\Delta_L \ge$ the bound below for every even $L \ge 4$ rungs. At $J_\perp = 2$ and $3/2$, every seed with $n_0 \in \{2, 4\}$ is rigorously rejected for all $\beta > 0$. The certificate is stored in [`results/certificate.json`](results/certificate.json), the route analysis is in [`IDEAS.md`](IDEAS.md), and the verification pipeline and failed perturbation route are documented in [`METHOD.md`](METHOD.md).

| $J_\perp$ | seed $(n_0, \beta_0)$ | certified gap $\ge$ | lengths covered |
|---|---|---:|---|
| 4 | $(4, 2)$ | 1.102403091 | every even $L \ge 4$ rungs ($N = 8, 12, 16, \dots$ spins) |
| 7/2 | $(4, 9/4)$ | 0.778648190 | every even $L \ge 4$ rungs |
| 3 | $(4, 5/2)$ | 0.454189334 | every even $L \ge 4$ rungs |
| 5/2 | $(4, 3)$ | 0.042433721 | every even $L \ge 4$ rungs |
| 163069/65536 | $(4, 3)$ | 0.001531165 | every even $L \ge 4$ rungs |
| 2, 3/2 | none ($n_0 \in \{2, 4\}$ rejected for all $\beta$) | — | no certificate |

## Model

The periodic two-leg ladder Hamiltonian on $L$ rungs is

$$H_L = \sum_{j=1}^{L} \sum_{a=1,2} \mathbf{S}_{a,j} \cdot \mathbf{S}_{a,j+1} + J_\perp \sum_{j=1}^{L} \mathbf{S}_{1,j} \cdot \mathbf{S}_{2,j},$$

with unit leg coupling and one rung ($d = 4$) per unit cell, so $L$ rungs contain $N = 2L$ spins. At $L = 2$, both periodic leg bonds are included (following the repository's two-cell ring convention). Because seed lengths $n_0$ are counted in rungs, a certificate at $n_0 = 4$ covers every even $L \ge 4$ rungs, corresponding to **$N = 8, 12, 16, \dots$ spins** ($N$ divisible by 4).

## Reflection positivity

In the standard spin-1/2 basis with real matrices $X = S^x$, $A = iS^y$, and $Z = S^z$, applying the checkerboard rotation $u = \exp(i\pi S^y) = 2A$ to one spin of every rung (alternating legs from rung to rung, which requires $L$ to be even) maps the Hamiltonian to a reflection-positive cell model with

$$-h_{\text{cell}} = \sum_\alpha C_\alpha \otimes C_\alpha, \qquad C_\alpha \in \{X \otimes I, A \otimes I, Z \otimes I, I \otimes X, I \otimes A, I \otimes Z\},$$

and real symmetric on-site rung term

$$o = -J_\perp (X \otimes X + A \otimes A + Z \otimes Z).$$

The function `exact_rp` in [`src/rp.jl`](src/rp.jl) verifies the rotation identities for both sublattice orientations over `Rational{BigInt}`, and `rp_model` stores the decomposition factors for verification by `CertifiedGaps.verify_rp` (in both `check.jl` and `check.jl --full`). Because the checkerboard rotation is a unitary equivalence, spectra are computed in the original $\mathrm{SU}(2)$ basis.

## Seeds and certified bounds

Each certified row establishes a unique ground state and $\Delta_L \ge$ the bound below for all even $L \ge 4$ rungs. The defects $\varepsilon_{\mathrm{s}}$ and $\varepsilon_\tau$ are exact rational upper bounds stored in the certificate (their decimal representations below are exact).

| $J_\perp$ | seed $(n_0, \beta_0)$ | $\varepsilon_{\mathrm{s}}$ | $\varepsilon_\tau$ | certified $\Delta \ge$ | on $\|J_\perp - J_*\| \le 10^{-4}$ |
|---|---|---:|---:|---:|---:|
| 4 | $(4, 2)$ | 0.03508321 | 0.02462420 | 1.102403091 | 1.063075969 |
| 7/2 | $(4, 9/4)$ | 0.05806853 | 0.02987064 | 0.778648190 | 0.751905791 |
| 3 | $(4, 5/2)$ | 0.09896358 | 0.04497730 | 0.454189334 | 0.435856503 |
| 5/2 | $(4, 3)$ | 0.18558434 | 0.05114736 | 0.042433721 | — |
| 163069/65536 $\approx 2.488235$ | $(4, 3)$ | 0.18772796 | 0.05269516 | 0.001531165 | — |
| 326133/131072 $\approx 2.488197$ | all 128 grid seeds rejected | | | — | |
| 2 | $n_0 \in \{2, 4\}$ rejected for every $\beta > 0$ | | | — | |
| 3/2 | $n_0 \in \{2, 4\}$ rejected for every $\beta > 0$ | | | — | |

As exact rationals, the $J_\perp = 4$ seed is $(n_0, \beta_0, \varepsilon_{\mathrm{s}}, \varepsilon_\tau) = (4, 2, 3508321/100000000, 123121/5000000)$. Asymmetric defect pairs where $\varepsilon_{\mathrm{s}}$ exceeds the symmetric diagonal threshold $\varepsilon_* \approx 0.12731$ (at $J_\perp \le 5/2$) are certified directly by the asymmetric two-purity recurrence in `certified_gap_bound`.

### Verification method

For each coupling, the suite encloses every eigenvalue of the $2$-, $4$-, and $8$-rung rings (using symmetry blocks with exact multiplicities, interval block entries, residual and orthogonality bounds, and Weyl's inequality), forms outward-rounded $\log Z$ sums, and classifies all 128 grid seeds $n_0 \in \{2, 4\}$, $\beta_0 = k/4$ ($k = 1, \dots, 64$) into three categories:
- `certified`: upper defect enclosures close the doubling recurrence;
- `rejected_by_recurrence`: lower defect bounds carry an exact rational witness showing that the recurrence diverges above $1/2$;
- `unresolved`: the enclosure is too wide to decide (does not occur for any grid point).

At $J_\perp = 2$ and $3/2$, monotone temperature-window bounds and a ground-state energy tail extend the rejection to every $\beta > 0$. Certified seeds at $J_\perp = 4, 7/2, 3$ are continued across $|J_\perp - J_*| \le 10^{-4}$ using an operator-norm bound on the centered rung sum (see [`METHOD.md`](METHOD.md)).

### Finite-grid threshold and negative rows

On the grid $n_0 \in \{2, 4\}$, $\beta_0 = k/4$ ($k = 1, \dots, 64$), every seed is rigorously rejected at $J_\perp = 326133/131072$ and a seed is certified at $J_\perp = 163069/65536$, locating the grid transition within $(2.48819732666015625, 2.4882354736328125)$. This brackets the threshold of the $(n_0 \le 4)$ seed grid rather than the physical gap-closing coupling (the Heisenberg ladder is gapped for all $J_\perp > 0$). Similarly, the all-$\beta$ rejections at $J_\perp = 2$ and $3/2$ show only that $n_0 \in \{2, 4\}$ is too small at those couplings; the heuristic fit in `results/scans/extrapolation.json` estimates that $n_0 = 6$ ($J_\perp = 2$) and $n_0 = 8$ ($J_\perp = 3/2$) would be required, corresponding to $24$- and $32$-spin rings.

## Reproduce

From the repository root, with one BLAS thread:

```sh
export OPENBLAS_NUM_THREADS=1
julia --project=models/heisenberg_ladder -e 'using Pkg; Pkg.instantiate()'
julia --project=models/heisenberg_ladder models/heisenberg_ladder/certify.jl
julia --project=models/heisenberg_ladder models/heisenberg_ladder/check.jl
julia --project=models/heisenberg_ladder models/heisenberg_ladder/check.jl --full
julia --project=models/heisenberg_ladder models/heisenberg_ladder/test/runtests.jl
julia --project=models/heisenberg_ladder models/heisenberg_ladder/scripts/audit.jl
julia --project=models/heisenberg_ladder models/heisenberg_ladder/scripts/replay_scans.jl
```

- `check.jl` (fast mode, 2 s): verifies the model-local source hashes, exact RP identities and the
  stored model, all grid decisions from stored rational defects (recurrence bounds and rejection
  witnesses), the best seeds, coupling continuations, and README gap claims (`CLAIMED_GAPS`)
  without recomputing spectra.
- `check.jl --full` (90 s): recomputes all verified spectra, defect and $\log Z$ enclosures, all
  1,024 grid decisions, and the all-$\beta$ exclusions from scratch.
- `certify.jl` (87.5 s on an idle machine, 10–30 s per coupling): generates
  `results/certificate.json`.
- `scripts/audit.jl` (11 s) and `scripts/replay_scans.jl` (1 s): the independent MPFR seed audit
  (see "What is rigorous") and the replay of the stored scans.
- `test/runtests.jl`: 6 s.

Run times are for Julia 1.12.4 on an Apple M5 with one BLAS thread, excluding package installation
and Julia start-up (wall-clock including start-up: `certify.jl` 93 s, or 175 s when rerun alongside
another job; `check.jl` 4 s; `check.jl --full` 93 s; `scripts/audit.jl` 17 s;
`scripts/replay_scans.jl` 4 s). An earlier run of the same suite on a slower Linux machine took
367 s to certify and 308 s to check from scratch.

The exploratory floating-point scans are regenerated per coupling by `scripts/scan.jl <J>` (20–35
min each, dominated by the 10-rung ring), followed by `scripts/extrapolate.jl`;
`scripts/threshold_bisection.jl` proposes the threshold endpoints. Running the six scans
concurrently once exceeded 17 GiB of RAM, so run at most one or two at a time.

## What is rigorous

- **Verified by `check.jl --full` from fresh computation:** exact RP identities, exact symmetry
  sector accounting, all eigenvalue and $\log Z$ enclosures, all 1,024 grid decisions with their
  upper/lower defect bounds and rejection witnesses, the core recurrence gap bounds, all-$\beta$
  exclusions, coupling continuations, README claims (`CLAIMED_GAPS`), and the SHA-256 hashes of the
  model-local sources `src/`, `certify.jl`, `check.jl`, and the model environment (tests and
  `scripts/` are not hashed; the core's hashes are recorded as provenance only, since the checker
  reruns the core). No stored spectrum or scan is an input. The fast `check.jl` re-establishes the
  parts that need no spectrum and trusts the stored defect enclosures.
- **Independent audit:** `scripts/audit.jl` redoes the $J_\perp = 4$ seed with 256-bit
  directed-rounding MPFR partition sums and checks the exact trace moments $\operatorname{Tr} H = 0$
  and $\operatorname{Tr}(H^2)$ against the eigenvalue enclosures.
- **Trusted base:** the two-purity doubling criterion (proved in `notes/` and formalized in `lean/`
  with the seeds as hypotheses) and IEEE binary64 round-to-nearest arithmetic in OpenBLAS/LAPACK
  under the explicit rounding-error allowances in [`METHOD.md`](METHOD.md).
- **Floating-point only (labeled as such, never used in proofs):** everything under `results/scans/`
  (dense scans up to 10 rungs / 20 spins / $2^{20}$ states, their $\log Z$ grids and finite-ring
  gaps, the threshold bisection, and the seed-scale extrapolation) and all timings.

## Layout

```text
certify.jl, check.jl             certificate generator and verifier (fast or --full)
METHOD.md                        certificate pipeline, perturbation analysis, and continuation
src/HeisenbergLadder.jl          module entry point (rp.jl, spectra.jl, suite.jl)
test/runtests.jl                 fast unit and verification tests
scripts/audit.jl                 independent MPFR audit of the J⊥ = 4 seed
scripts/                         exploratory floating-point scans and replay (not certificates)
results/certificate.json         stored certificate (single-line JSON, exact rational strings)
results/audit.json               output of scripts/audit.jl
results/scan_replay.json         output of scripts/replay_scans.jl
results/scans/                   floating-point scan data
```

Hash pinning: `results/certificate.json` records by repository path the SHA-256 hashes of the
model-local sources `src/`, `certify.jl`, `check.jl`, `Project.toml`, and `Manifest.toml`
(`metadata.source_sha256`), and both checker modes require exactly this key set and every digest.
The shared core (`julia/src`, `julia/Project.toml`, `julia/Manifest.toml`) is recorded under
`metadata.core_sha256` as provenance only and is not compared with the checkout; `notes/` and
`METHOD.md` are not hashed. The certificate is regenerated by `certify.jl` in under two minutes, so
there are no layout exceptions.

## History

During the directory-format standardization, `PERTURBATION.md` was renamed to `METHOD.md`
(incorporating the pipeline documentation from this README), `audit.jl` was moved to `scripts/`,
`check.jl` was split into a fast default mode and `--full`, and `results/certificate.json` and
`results/audit.json` were regenerated because `check.jl` and `src/suite.jl` are hashed. Every
certified seed, defect, decision, and gap bound is identical to the previous certificate.

The source-hash policy change (the checkers enforce only model-local sources and their key set;
core hashes are provenance; the scope note is no longer hashed) regenerated
`results/certificate.json`, `results/audit.json`, and `results/scan_replay.json`. A deep comparison
of old and new JSON found every enclosure, defect, decision, exclusion, continuation, and bound
identical; only hashes, timestamps, and timings changed.
