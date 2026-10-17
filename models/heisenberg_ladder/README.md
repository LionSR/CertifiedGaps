# Two-leg spin-1/2 Heisenberg ladder

This module certifies uniform spectral gaps for the two-leg spin-1/2 Heisenberg ladder at rung couplings $J_\perp = 4$, $7/2$, $3$, $5/2$, and $163069/65536$ from verified complete spectra of rings of $2$, $4$, and $8$ rungs, computed entirely with the shared core pipeline (`CertifiedGaps.GapCertificates.certify_gap`). Each certificate establishes a unique ground state and $\Delta_L \ge$ the bound below for every even $L \ge 4$ rungs. At $J_\perp = 2$ and $3/2$, every seed with $n_0 \in \{2, 4\}$ is rigorously rejected for all $\beta > 0$. The certificates are `results/certificate_J*.json` (one core certificate per coupling), the grid decisions and negative results are in [`results/suite.json`](results/suite.json), the route analysis is in [`IDEAS.md`](IDEAS.md), and the method and the failed perturbation route are documented in [`METHOD.md`](METHOD.md).

| $J_\perp$ | seed $(n_0, \beta_0)$ | certified gap $\ge$ | lengths covered |
|---|---|---:|---|
| 4 | $(4, 2)$ | 1.102404331 | every even $L \ge 4$ rungs ($N = 8, 12, 16, \dots$ spins) |
| 7/2 | $(4, 9/4)$ | 0.778648930 | every even $L \ge 4$ rungs |
| 3 | $(4, 5/2)$ | 0.454189766 | every even $L \ge 4$ rungs |
| 5/2 | $(4, 3)$ | 0.042434525 | every even $L \ge 4$ rungs |
| 163069/65536 | $(4, 3)$ | 0.001540040 | every even $L \ge 4$ rungs |
| 2, 3/2 | none ($n_0 \in \{2, 4\}$ rejected for all $\beta$) | — | no certificate |

The gap column is the core interval bound `interval_bound_decimal` (`certified_gap_bound`, rounded down to nine digits); the exact rational chain bound of each certificate is listed below.

## Model

The periodic two-leg ladder Hamiltonian on $L$ rungs is

$$H_L = \sum_{j=1}^{L} \sum_{a=1,2} \mathbf{S}_{a,j} \cdot \mathbf{S}_{a,j+1} + J_\perp \sum_{j=1}^{L} \mathbf{S}_{1,j} \cdot \mathbf{S}_{2,j},$$

with unit leg coupling and one rung ($d = 4$) per unit cell, so $L$ rungs contain $N = 2L$ spins. At $L = 2$, both periodic leg bonds are included (following the repository's two-cell ring convention). Because seed lengths $n_0$ are counted in rungs, a certificate at $n_0 = 4$ covers every even $L \ge 4$ rungs, corresponding to **$N = 8, 12, 16, \dots$ spins** ($N$ divisible by 4).

## Reflection positivity

Write $X = S^x$, $A = iS^y$ (real antisymmetric), and $Z = S^z$, so $\mathbf{S} \cdot \mathbf{S} = X \otimes X - A \otimes A + Z \otimes Z$. The certificate is stated in the **leg-rotated frame**
$U_L = \bigotimes_{j} (1 \otimes u)$, $u = \exp(i\pi S^y) = 2A$, which rotates the leg-2 spin of every rung (uniformly, for every $L$). The leg bonds are unchanged, the rung term becomes

$$o = U (J_\perp \mathbf{S}_1 \cdot \mathbf{S}_2) U^\mathsf{T} = -J_\perp (X \otimes X + A \otimes A + Z \otimes Z),$$

and the physical $S^z_{\text{tot}}$ becomes the cell charge $z_1 - z_2$. With the internal cell reflection $F = \sigma^y \otimes \sigma^y = -4\,A \otimes A$ (a real symmetric involution with $F o F = o$), the bond is reflection positive in the reflected form

$$-h = \sum_\alpha C_\alpha \otimes F \overline{C_\alpha} F, \qquad C_\alpha \in \{X \otimes 1, A \otimes 1, Z \otimes 1, 1 \otimes X, 1 \otimes A, 1 \otimes Z\},$$

since $\sigma^y X \sigma^y = -X$, $\sigma^y A \sigma^y = A$, $\sigma^y Z \sigma^y = -Z$ on each leg. The core `RingModel` (`ladder_model` in [`src/model.jl`](src/model.jl)) carries the charge $z_1 - z_2$, the charge flip $F$, and the on-site $\mathbb{Z}_2$ symmetry leg exchange $U\,\mathrm{SWAP}\,U^\mathsf{T}$; the core checks all of them exactly. `ladder_rp` stores the factors and $F$ for `CertifiedGaps.verify_rp`, and `check_frame` verifies over `Rational{BigInt}` that $U$ maps the physical bond, rung term, and $S^z_{\text{tot}}$ to the rotated ones (both checker modes and `test/runtests.jl`, which also compares the spectra of both frames with an independent dense Hamiltonian).

## Seeds and certified bounds

Each certified coupling has one core certificate `results/certificate_J<tag>.json` with the exact seed $(n_0, \beta_0, \varepsilon_{\mathrm{s}}, \varepsilon_\tau)$ (defects rounded up to $10^{-15}$; decimals below are exact), the exact rational chain, and both gap bounds.

| $J_\perp$ | seed $(n_0, \beta_0)$ | $\varepsilon_{\mathrm{s}}$ | $\varepsilon_\tau$ | interval $\Delta \ge$ | exact chain $\Delta >$ | on $\|J_\perp - J_*\| \le 10^{-4}$ |
|---|---|---:|---:|---:|---:|---:|
| 4 | $(4, 2)$ | 0.035083179771961 | 0.024624083715215 | 1.102404331 | 1.101612406 | 1.063077260 |
| 7/2 | $(4, 9/4)$ | 0.058068500793985 | 0.029870522063335 | 0.778648930 | 0.777944996 | 0.751906573 |
| 3 | $(4, 5/2)$ | 0.098963550759311 | 0.044977192016543 | 0.454189766 | 0.453872996 | 0.435856990 |
| 5/2 | $(4, 3)$ | 0.185584310567304 | 0.051147246075377 | 0.042434525 | 0.042401527 | — |
| 163069/65536 $\approx 2.488235$ | $(4, 3)$ | 0.187727933251953 | 0.052695044875928 | 0.001540040 | 0.001535915 | — |
| 326133/131072 $\approx 2.488197$ | all 128 grid seeds rejected | | | — | | |
| 2 | $n_0 \in \{2, 4\}$ rejected for every $\beta > 0$ | | | — | | |
| 3/2 | $n_0 \in \{2, 4\}$ rejected for every $\beta > 0$ | | | — | | |

Asymmetric defect pairs where $\varepsilon_{\mathrm{s}}$ exceeds the symmetric diagonal threshold $\varepsilon_* \approx 0.12731$ (at $J_\perp \le 5/2$) are certified directly by the asymmetric two-purity recurrence.

### Verification method

For each coupling, `verified_ring_spectrum` encloses every eigenvalue of the $2$-, $4$-, and $8$-rung rings in blocks of charge, momentum, leg-exchange parity, and (at $S^z = 0$) flip parity, with exact multiplicities (largest block 715 at 8 rungs). `seed_bounds` turns them into defect enclosures for all 128 grid seeds $n_0 \in \{2, 4\}$, $\beta_0 = k/4$ ($k = 1, \dots, 64$), reusing the spectra, and each seed is classified as
- `certified`: the upper defects close the recurrence (`certified_gap_bound`);
- `rejected_by_recurrence`: the lower defects carry an exact rational witness that the recurrence diverges above $1/2$;
- `unresolved`: neither (does not occur).

The best certified seed is passed to `certify_gap`, which runs the exact rational chain. At $J_\perp = 2$ and $3/2$, monotone temperature-window bounds and a ground-energy tail (from `shifted_log_partition` and `ground_energy_enclosure`) extend the rejection to every $\beta > 0$. Certified seeds at $J_\perp = 4, 7/2, 3$ are continued across $|J_\perp - J_*| \le 10^{-4}$ using an operator-norm bound on the centered rung sum (see [`METHOD.md`](METHOD.md)).

### Finite-grid threshold and negative rows

On the grid $n_0 \in \{2, 4\}$, $\beta_0 = k/4$ ($k = 1, \dots, 64$), every seed is rigorously rejected at $J_\perp = 326133/131072$ and a seed is certified at $J_\perp = 163069/65536$, locating the grid transition within $(2.48819732666015625, 2.4882354736328125)$. This brackets the threshold of the $(n_0 \le 4)$ seed grid rather than the physical gap-closing coupling (the Heisenberg ladder is gapped for all $J_\perp > 0$). Similarly, the all-$\beta$ rejections at $J_\perp = 2$ and $3/2$ show only that $n_0 \in \{2, 4\}$ is too small at those couplings; the heuristic fit in `results/scans/extrapolation.json` estimates that $n_0 = 6$ ($J_\perp = 2$) and $n_0 = 8$ ($J_\perp = 3/2$) would be required, corresponding to $24$- and $32$-spin rings.

## Reproduce

From the repository root, in the core environment with one BLAS thread:

```sh
export OPENBLAS_NUM_THREADS=1
julia --project=julia -e 'using Pkg; Pkg.instantiate()'
julia --project=julia models/heisenberg_ladder/certify.jl
julia --project=julia models/heisenberg_ladder/check.jl
julia --project=julia models/heisenberg_ladder/check.jl --full
julia --project=julia models/heisenberg_ladder/test/runtests.jl
julia --project=julia models/heisenberg_ladder/scripts/audit.jl
julia --project=julia models/heisenberg_ladder/scripts/replay_scans.jl
```

- `check.jl` (fast): model-local source hashes (`Results.verify_model_metadata`), exact frame and
  RP identities, all grid decisions from the stored rational defects, the best seeds, each core
  certificate with `check_certificate(cert; spectra = false)` (RP and the exact chain), the
  continuations, and the README claims (`CLAIMED_GAPS`), without spectra.
- `check.jl --full`: also recomputes every verified spectrum, all 1,024 grid defect enclosures,
  each certificate's seed (`check_certificate` with the fresh spectra), and the all-$\beta$
  exclusions.
- `certify.jl`: writes the five core certificates and `results/suite.json`.
- `scripts/audit.jl` and `scripts/replay_scans.jl`: the independent MPFR audit of the
  $J_\perp = 4$ seed and the replay of the stored scans.

Measured CPU (user) times, Julia 1.12 on an Apple M5 with one BLAS thread, including start-up:
`certify.jl` 120 s (156 s wall clock on a lightly loaded machine; 15–20 s per coupling, about 13 s
of it the 8-rung spectrum), `check.jl` 19 s, `check.jl --full` 116 s, `test/runtests.jl` 30 s,
`scripts/audit.jl` 30 s, `scripts/replay_scans.jl` 9 s. The regeneration ran alongside many other
jobs (load average 70–130), so wall-clock times were up to 3.6 times larger.

The exploratory floating-point scans are regenerated per coupling by `scripts/scan.jl <J>`
(unverified eigenvalues of the core symmetry blocks up to 10 rungs), followed by
`scripts/extrapolate.jl`; `scripts/threshold_bisection.jl` proposes the threshold endpoints. The
stored scans were produced by the earlier local block builder and were not rerun.

## What is rigorous

- **Verified by `check.jl --full` from fresh computation:** exact frame and RP identities, the
  exact symmetry blocks and multiplicities (`ring_blocks`), all eigenvalue and $\log Z$
  enclosures, all 1,024 grid decisions with their upper/lower defect bounds and rejection
  witnesses, each core certificate's seed and exact rational chain, the core interval bounds, the
  all-$\beta$ exclusions, the coupling continuations, the README claims (`CLAIMED_GAPS`), and the
  SHA-256 hashes of the model-local sources `src/`, `certify.jl`, and `check.jl` (tests and
  `scripts/` are not hashed; the core and its environment are recorded as provenance only, since
  the checker reruns the core). No stored spectrum or scan is an input. The fast `check.jl`
  re-establishes the parts that need no spectrum and trusts the stored defect enclosures.
- **Enclosure:** each block's Float64 center (exactly Hermitian, real when the block is real) and an
  exact Frobenius bound on its distance to the interval block go to
  `CertifiedGaps.VerifiedSpectra.verified_eigen`: residual and orthogonality bounds evaluated in
  Float64 (BLAS allowed) with the a priori allowances $\gamma_{N+8}$ (real) or $\gamma_{2N+8}$
  (complex, four real products) on Frobenius norms, then the smaller of the Ostrowski and
  polar-factor radii plus Weyl for the block error, and outward-rounded $\log Z$ sums
  (`log_partition_enclosure`).
- **Independent audit:** `scripts/audit.jl` redoes the $J_\perp = 4$ seed with 256-bit
  directed-rounding MPFR partition sums and checks the exact trace moments $\operatorname{Tr} H = 0$
  and $\operatorname{Tr}(H^2)$ against the eigenvalue enclosures.
- **Trusted base:** the two-purity doubling criterion with an internal cell reflection (proved in
  `notes/` and formalized in `lean/` with the seeds as hypotheses) and IEEE binary64
  round-to-nearest arithmetic with gradual underflow in OpenBLAS/LAPACK under the $\gamma_n$
  rounding model, valid for any summation order and FMA.
- **Floating-point only (labeled as such, never used in proofs):** everything under `results/scans/`
  (dense scans up to 10 rungs / 20 spins / $2^{20}$ states, their $\log Z$ grids and finite-ring
  gaps, the threshold bisection, and the seed-scale extrapolation) and all timings.

## Layout

```text
certify.jl, check.jl             certificate generator and verifier (fast or --full)
METHOD.md                        frame, pipeline, perturbation analysis, and continuation
src/HeisenbergLadder.jl          module entry point (model.jl, suite.jl)
src/model.jl                     core RingModel in the leg-rotated frame, RP data, frame identities
src/suite.jl                     grid, rejections, exclusions, continuation, produce and check
test/runtests.jl                 fast unit and verification tests
scripts/audit.jl                 independent MPFR audit of the J⊥ = 4 seed
scripts/                         exploratory floating-point scans and replay (not certificates)
results/certificate_J*.json      core gap certificates (J4, J7_2, J3, J5_2, J163069_65536)
results/suite.json               grid decisions, rejection witnesses, exclusions, continuations
results/audit.json               output of scripts/audit.jl
results/scan_replay.json         output of scripts/replay_scans.jl
results/scans/                   floating-point scan data
```

Hash pinning: every result records by repository path the SHA-256 hashes of `src/`, `certify.jl`,
and `check.jl` (`metadata.source_sha256`), and both checker modes require exactly this key set and
every digest. The shared core and its environment (`julia/src`, `julia/Project.toml`,
`julia/Manifest.toml`) are recorded under `metadata.core_sha256` as provenance only; `notes/` and
`METHOD.md` are not hashed. The certificates are regenerated by `certify.jl` in about two
minutes, so there are no layout exceptions.

## History

During the directory-format standardization, `PERTURBATION.md` was renamed to `METHOD.md`,
`audit.jl` moved to `scripts/`, and `check.jl` was split into a fast default mode and `--full`;
later, the source-hash policy change and the move of duplicated helpers into the core regenerated
the results with every certified number identical.

The certificate pipeline then moved entirely into the core: the local block builder, eigenvalue
enclosure ($\gamma_{8n+32}$ with $\|H\|_F$), $\log Z$ sums, and the single-file
`results/certificate.json` were replaced by `verified_ring_spectrum`, `seed_bounds`, and
`certify_gap` in the uniform leg-rotated frame (cell reflection $\sigma^y \otimes \sigma^y$,
conserved $S^z$, flip and leg-exchange parities), with one core certificate per coupling and the
model environment replaced by the core's. Every grid decision, best seed, and exclusion is
unchanged, every new defect bound is at most the old one, and the interval gap bounds rose
(old → new): $1.102403091 \to 1.102404331$, $0.778648190 \to 0.778648930$,
$0.454189334 \to 0.454189766$, $0.042433721 \to 0.042434525$, $0.001531165 \to 0.001540040$.
