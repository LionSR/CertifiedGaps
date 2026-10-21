# Spin-1 large-$`D`$ chain: uniform-gap certificates

We certify uniform spectral gaps for the periodic spin-1 Heisenberg chain with single-ion
anisotropy $D(S^z)^2$ in two stages. Stage 1 runs the shared core pipeline
(`CertifiedGaps.GapCertificates.certify_point`) on complete verified spectra of the rings $L = 6$
and $L = 12$. Stage 2 reuses those stored spectra for spatial transfer-operator squaring, an exact
purity coupling identity that eliminates the 24-site ring, and centered Weyl perturbation transport
in $D$. Every certificate gives a unique ground state and a uniform gap
$\Delta_L \ge \Delta_{\mathrm{cert}}$ for every even $L \ge 6$ sites. Decimal bounds are rounded
**down**; exact rationals are stored in `results/`.

| Single-ion anisotropy $D$ | Seed $(n_0, \beta_0)$ | Certified gap $\Delta_L \ge$ | Lengths covered |
|---|---|---:|---|
| $4$ | $(12, 29/10)$ + finite bridge at $L = 6, 8, 10$ | **$0.899422250$** | every even $L \ge 6$ sites |
| $[3.2059998, 3.206]$ | $(6, 381/125)$, Weyl-transported | **$0.001823729$** | every even $L \ge 6$ sites |

## Model

```math
H_L = \sum_{j=1}^L \mathbf{S}_j \cdot \mathbf{S}_{j+1} + D \sum_{j=1}^L (S_j^z)^2, \qquad \mathbf{S}_j^2 = 2I.
```
The unit cell is one spin-1 site ($c = 1$, local dimension $3$), so all lengths $L$ and $n_0$
count sites. The core `RingModel` (`ring_model` in `src/model.jl`) uses the $S^z$ basis
$(|{+1}\rangle, |0\rangle, |{-1}\rangle)$ with the integer matrices of
`CertifiedGaps.Spins.integer_spin_one`, $`x = \sqrt{2}\,S^x`$ and the real antisymmetric
$`a = \sqrt{2}\,iS^y`$:
```math
h = S^z \otimes S^z + \tfrac{1}{2}\, x \otimes x - \tfrac{1}{2}\, a \otimes a, \qquad o = D\mathrm{diag}(1, 0, 1),
```
with conserved charge $`m \in \{1, 0, -1\}`$, the spin flip $F = \mathrm{antidiag}(1, 1, 1)$, and
the lattice reflection $j \mapsto -j$ with trivial on-site map ($P = 1$, since $h$ is symmetric
under the swap of its two sites).

**Frame.** The Lean model (`code/lean/CertifiedGaps/Models/SpinOneLargeD.lean`) uses the Cartesian
basis with $`(S^\alpha)_{ab} = -i\varepsilon_{\alpha ab}`$. The unitary $U = W N^{-1/2}$ with
```math
W = \begin{pmatrix} -i & 0 & i \\ 1 & 0 & 1 \\ 0 & i & 0 \end{pmatrix}, \qquad N = W^\dagger W = \mathrm{diag}(2, 1, 2)
```
maps it to the $S^z$ basis: $U^\dagger S^x U = x/\sqrt{2}$, $U^\dagger (iS^y) U = a/\sqrt{2}$,
$U^\dagger S^z U = S^z$ and $U^\dagger (S^z)^2 U = \mathrm{diag}(1, 0, 1)$, so $U^{\otimes L}$
conjugates one ring Hamiltonian into the other. The $\sqrt{2}$ in $U$ rules out the core's exact
`is_frame_change`, so `cartesian_frame` checks the equivalent identities for $W$
($W^\dagger S^x W = x$, $W^\dagger (iS^y) W = a$, $W^\dagger S^z W = N S^z$,
$W^\dagger (S^z)^2 W = N (S^z)^2$) in exact Gaussian-integer arithmetic, and `test/frame.jl`
also compares the verified spectra with dense diagonalization of the Cartesian Hamiltonian for
$L \le 6$.

## Reflection positivity

With the real symmetric involution $R = \mathrm{antidiag}(1, -1, 1)$ (`RP_REFLECTION`), which
commutes with $o$, one has $R S^z R = -S^z$, $R x R = -x$ and $R a R = a$, hence
```math
-h = S^z \otimes R S^z R + \tfrac{1}{2}\, x \otimes R x R + \tfrac{1}{2}\, a \otimes R a R,
```
the "reflected" RP form with integer factors and weights $(1, 1/2, 1/2)$ (`rp_data`). It is verified
exactly by `CertifiedGaps.RingModels.verify_rp` in `certify_gap`, in `check.jl` (`check_point`), and
in `test/frame.jl` (which also checks that the exact LDL factorization of the realignment finds an
RP form with $R$ and that the bond is not RP in this basis without $R$). In the Cartesian frame it
is the familiar $-h = \sum_\alpha G_\alpha \otimes G_\alpha$ with $G_\alpha = i S^\alpha$ real and
antisymmetric.

## Seeds and certified bounds

Stage-1 core certificates (`results/certificate_D*.json`, $n_0 = 6$). The table reports the
interval bound `interval_bound_decimal` (`certified_gap_bound`) and the exact rational chain
bound `gap_lower_bound`; both are rigorous.

| $D$ | $\beta_0$ | $\varepsilon_{\mathrm{s}}$ | $\varepsilon_\tau$ | gap $\ge$ (interval) | gap $>$ (chain) |
|---|---|---|---|---:|---:|
| $4$ | $9/4$ | `64102636591415//10^15` | `40520211376808//10^15` | $0.694746846$ | $0.694042913$ |
| $1603/500$ | $381/125$ | `191641808998499//10^15` | `46451711086584//10^15` | $0.002754398$ | $0.002750338$ |

Both seeds lie below the defects assumed by the Lean file ($\varepsilon_{\mathrm{s}} \le 0.06411$,
$\varepsilon_\tau \le 0.04053$ at $D = 4$; $\varepsilon_{\mathrm{s}} \le 0.19164415$,
$\varepsilon_\tau \le 0.046462453$ at $D = 1603/500$), which `check.jl` verifies exactly.

Stage-2 cases (`results/certificate.json`, interval bounds, every even $L \ge n_0$):

| Case | $D$ | $(n_0, \beta_0)$ | $\varepsilon_{\mathrm{s}}$ | $\varepsilon_\tau$ | gap $\ge$ |
|---|---|---|---:|---:|---:|
| reference | $4$ | $(6, 9/4)$ | $0.0641027$ | $0.0405203$ | $0.694746846$ |
| survey seed | $4$ | $(6, 3)$ | $0.0809098$ | $0.00643407$ | $0.543857041$ |
| space lift 1 | $4$ | $(12, 29/10)$ | $0.00363252$ | $0.0314054$ | **$0.899422250$** |
| space lift 2 | $4$ | $(24, 121/40)$ | $7.81249 \times 10^{-6}$ | $0.0565974$ | $0.660562073$ |
| space lift 3 | $4$ | $(48, 63/20)$ | $4.18929 \times 10^{-11}$ | $0.104277$ | $0.403276911$ |
| space lift 4 | $4$ | $(96, 133/40)$ | $2.07960 \times 10^{-21}$ | $0.192979$ | $0.064969369$ |
| transport | $[3.9979773, 4]$ | $(6, 23/10)$ | $0.115897$ | $0.137682$ | $0.002594199$ |
| reference | $1603/500$ | $(6, 381/125)$ | $0.191642$ | $0.0464518$ | $0.002754398$ |
| transport | $[3.2059998, 3.206]$ | $(6, 381/125)$ | $0.191648$ | $0.0464657$ | **$0.001823729$** |

At $D = 4$ the $n_0 = 12$ bound extends to every even $L \ge 6$ through the finite bridge: the
verified gaps of the rings $L = 6, 8, 10$ (`results/finite_bridge.json`) exceed it. For
comparison, $\Delta_{12} \in [2.2694582, 2.2694597]$ at $D = 4$ and
$\Delta_{12} \in [1.5453276, 1.5453288]$ at $D = 3.206$.

- **Stage 1 (`certify.jl`):** total $S^z$, momentum $0 \le m \le 6$, and flip parity at $S^z = 0$
  (the core `ring_blocks` with `flip` and `reflection`) split the $3^{12} = 531{,}441$ states into
  92 blocks of dimension at most $5{,}798$; the lattice reflection makes every block a real
  symmetric matrix. Each block is enclosed by `CertifiedGaps.VerifiedSpectra.verified_eigen`
  (largest radius $5.2 \times 10^{-7}$), and $\log Z$ is summed by `log_partition_enclosure`. The
  spectra are stored with `write_spectrum`, and their SHA-256 is pinned in the certificate
  (`spectrum_sha256`).
- **Stage 2 (`scripts/transfer_certify.jl`):** from the stored $L = 6, 12$ spectra at $\beta$ and
  $2\beta$, transfer-operator squaring bounds $P_{\mathrm{s}}(12, \beta_0)$, and the exact coupling
  identity
  ```math
  P_\tau(4n, \beta_0) = \frac{P_{\mathrm{s}}(2n, 2\beta_0)\,P_\tau(2n, \beta_0)^2}{P_{\mathrm{s}}(2n, \beta_0)^2}
  ```
  bounds $P_\tau(24, \beta_0)$ without the 24-site ring (`paper/thermal-mpo-gap-scope.tex`,
  Appendix E). A centered Weyl bound ($0 \le \sum_j (S_j^z)^2 \le L$) transports seeds to the
  fixed intervals $D \in [3.2059998, 3.206]$ and $[3.9979773, 4]$. As a diagnostic, the driver
  also records how far bisection on a $10^{-7}$ grid could extend these intervals with the
  current enclosures (to $3.2059997$ and $3.9979771$). Proofs are in [`METHOD.md`](METHOD.md),
  and the route evaluation is in [`IDEAS.md`](IDEAS.md).

### Exploratory numerical scan (floating point, non-certifying)

`results/scan_summary.json` and the numerical points of `results/transfer_curve.*` come from a
floating-point scan of $`D \in \{1, 1.25, \dots, 4\} \cup \{3.205, 3.206\}`$,
$`n_0 \in \{2, \dots, 6\}`$ and $`\beta_0 \in \{0.125, \dots, 20\}`$ (refined to $2.75:0.002:3.75$ at
$D = 3.205, 3.206$):

| $D$ | Best even seed $(n_0, \beta_0)$ | Numerical bound | $\Delta_{12}$ (ED) | $\max_\beta \min(P_{\mathrm{s}}, P_\tau)$ at $n_0 = 6$ | $\delta(6) = E_0(12) - 2E_0(6)$ |
|---:|---|---:|---:|---:|---:|
| $1.0$ | none at $n_0 \le 6$ | — | $0.2554$ | $0.335$ | $0.3528$ |
| $1.5$ | none at $n_0 \le 6$ | — | $0.3418$ | $0.375$ | $0.3229$ |
| $2.0$ | none at $n_0 \le 6$ | — | $0.5840$ | $0.494$ | $0.2422$ |
| $2.5$ | none at $n_0 \le 6$ | — | $0.9452$ | $0.655$ | $0.1505$ |
| $3.0$ | none at $n_0 \le 6$ | — | $1.3641$ | $0.799$ | $0.0857$ |
| $3.205$ | none at $n_0 \le 6$ | — | $1.5444$ | $0.839$ | $0.0677$ |
| $3.206$ | $(6, 3.048)$ | $0.00283$ | $1.5453$ | $0.839$ | $0.0676$ |
| $3.25$ | $(6, 3.0)$ | $0.0931$ | $1.5845$ | $0.847$ | $0.0642$ |
| $3.5$ | $(6, 2.625)$ | $0.3232$ | $1.8095$ | $0.888$ | $0.0483$ |
| $3.75$ | $(6, 2.375)$ | $0.5121$ | $2.0381$ | $0.916$ | $0.0365$ |
| $4.0$ | $(6, 2.25)$ | $0.6948$ | $2.2695$ | $0.938$ | $0.0278$ |

At $n_0 = 6$, $D = 3.206$ sits at the boundary of the admissible region for direct seeds; smaller
$D$ needs larger $n_0$. A single $L = 14$, $S^z = 1$ momentum block has dimension $41{,}835$
($14.0$ GB as one dense real matrix, about eight times that for a verified enclosure; see
`results/feasibility.json`).

## Reproduce

From the repository root, in the shared environment (`OPENBLAS_NUM_THREADS=1`):

```sh
export OPENBLAS_NUM_THREADS=1
julia --project=code/julia code/models/spin1_largeD/certify.jl              # Stage 1 (see below)
julia --project=code/julia code/models/spin1_largeD/scripts/transfer_certify.jl   # Stage 2 (~78 s)
julia --project=code/julia code/models/spin1_largeD/scripts/report.jl --certified-only   # curve points
julia --project=code/julia code/models/spin1_largeD/check.jl                # fast check (~44 s)
julia --project=code/julia code/models/spin1_largeD/check.jl --full         # + L = 6, 12 replay (~1 h 45 min)
julia --project=code/julia code/models/spin1_largeD/test/runtests.jl        # tests + fast check (~63 s)
```

Stage 1 costs about 1,800 s of CPU time per point (the $L = 12$ ring dominates), with a peak
resident memory of 2.2 GB. `certify.jl 1603//500` runs one point; the two points ran as separate
processes in 52 min wall time each on an Apple M-series host under a load average near 200, so
all times above are inflated. The full replay recomputes the same spectra as `certify.jl` (about
52 min per point under that load, run one point after the other); it was not rerun for this
revision. The floating-point scan uses the core's unverified block eigenvalues
(`float_ring_spectrum`):

```sh
julia --project=code/julia code/models/spin1_largeD/scripts/scan.jl [D ...]
julia --project=code/julia code/models/spin1_largeD/scripts/refine_beta.jl 1603//500 641//200
julia --project=code/julia code/models/spin1_largeD/scripts/report.jl
julia --project=code/julia code/models/spin1_largeD/scripts/feasibility.jl
```

## What is rigorous

- **Rigorous (relative to the trusted base):** the exact RP identity and frame change; the exact
  symmetry reduction and block completeness ($3^L$ states); every eigenvalue enclosure at
  $L = 6, 8, 10, 12$ by `CertifiedGaps.VerifiedSpectra.verified_eigen` under the IEEE binary64 +
  BLAS model (round to nearest, gradual underflow, the $\gamma_{N+8}$ dot-product bound for the
  real blocks, valid for any summation order and FMA, Frobenius norms, then Ostrowski or
  polar-factor Weyl bounds); outward-rounded $\log Z$ intervals over complete spectra; exact
  rational seed defects; transfer squaring, the purity coupling, and centered Weyl transport; and
  the recurrence bounds (`certified_gap_bound`, and the exact chain `rational_gap_chain` for
  stage 1).
- **Fast vs. full verification:** `check.jl` checks each core certificate with
  `GapCertificates.check_point` (stage-1 source hashes; exactly the point's model, RP data, seed
  and params; the Lean seed hypotheses read by `lean_points` from
  `code/lean/CertifiedGaps/Models/SpinOneLargeD.lean`; the exact chain and interval bound), then the
  SHA-256 of the stored spectra and the $\log Z$ and defects recomputed from them by
  `seed_bounds(...; spectra)`, a fresh $L = 6$ replay, every stage-2 case (exact equality), and a
  fresh $L = 8, 10$ bridge. `check.jl --full` re-verifies both points from the model
  (`check_point(...; spectra = true)`, dense $L = 6$ and $L = 12$).
- **Numerical only:** the $(D, n_0, \beta_0)$ scan, the $\Delta_{12}$ plot denominators of the
  scan curve, and the seed-size extrapolations for $D \le 3.2$.

## Layout

Stored results live in [`data/models/spin1_largeD/`](../../../data/models/spin1_largeD/);
this directory's `results/` is a relative symlink to that canonical location.

| Path | Role |
|---|---|
| `certify.jl` | Stage 1: core certificates and stored spectra at both points |
| `check.jl` | Verifier for both stages (`--full` re-verifies the spectra) |
| `src/model.jl` | Exact model, frame, RP data, points, Lean file path, file names |
| `src/TransferCertificate.jl` | Stage 2: transfer squaring, coupling, Weyl transport, bridge, checker |
| `scripts/transfer_certify.jl` | Stage 2 driver: `results/certificate.json`, `results/finite_bridge.json` |
| `scripts/` (others) | Floating-point scan, $\beta$ refinement, $L = 14$ feasibility, report/plot |
| `test/` | Frame and RP tests, stage-2 inequalities and mutation tests, then `check.jl` |
| `results/certificate_D*.json` | Stage-1 core certificates |
| `results/spectrum_D*_L{6,12}.json.gz` | Verified spectra (`write_spectrum`), pinned by SHA-256 |
| `results/` (others) | Stage-2 certificate, finite bridge, scan summary, curve, feasibility |

Provenance: stage-1 certificates record `model_metadata` of `src/model.jl` and `certify.jl`;
stage-2 files record `src/`, `scripts/transfer_certify.jl` and `check.jl`. The checker enforces
these model-local hashes; the core hashes are recorded as provenance only. No file is hash-pinned
by path beyond these sources and the stored spectra.

## History

**Core drivers and the lattice reflection.** `certify.jl` and `check.jl` now use the core
`certify_point`, `check_point` and `lean_points` (replacing hand-written Lean seed constants),
the spin matrices come from `Spins.integer_spin_one`, stage 2 bounds its defects with
`Seeds.seed_defects`, and the scan uses `float_ring_spectrum`, `seed_record` and `purity_peaks`.
The model gained the lattice reflection ($P = 1$), so every $L = 12$ block is real (same 92
blocks), and the largest radius fell from $1.54 \times 10^{-6}$ to $5.2 \times 10^{-7}$. All
certificates were regenerated with the same seeds, $\beta$ choices and intervals; every defect
decreased and every bound increased (old $\to$ new):

- stage 1 (interval): $0.694746580 \to 0.694746846$ at $D = 4$ and
  $0.002752708 \to 0.002754398$ at $D = 3.206$;
- the $n_0 = 12$ seed: $0.899422073 \to 0.899422250$;
- the transported intervals: $0.001821437 \to 0.001823729$ on $[3.2059998, 3.206]$ and
  $0.002590164 \to 0.002594199$ on $[3.9979773, 4]$.

The finite bridge is unchanged on its $10^{-6}$ grid.

**Consolidation onto the shared core.** Before that, the model had its own pipeline: a local block
builder (`model.jl`) with spin inversion, a local residual/polar enclosure with a priori
$\gamma_{8m+32}$ constants and $`\|H\|_2 \le L(|D| + 2)`$, run-length-encoded $L = 12$ spectra
(`src/spectra.jl`), a local $\log Z$ with a truncated tail, a per-model environment, and hash-pinned
top-level scripts. It was replaced by the shared core. The model is now a `RingModel` in the $S^z$
frame (with the Cartesian frame checked exactly). Stage 1 is `certify_gap` with `flip`, which gives
the same 92 blocks; the spectra are stored with `write_spectrum`; and stage 2 sums $\log Z$ over the
complete stored spectra. Both points were regenerated with the same seeds and $\beta$ choices. Every
defect decreased and every bound increased (old $\to$ new):

- stage 1: $0.694694192 \to 0.694746580$ at $D = 4$ and $0.002213456 \to 0.002752708$ at
  $D = 3.206$;
- the $n_0 = 12$ seed: $0.899143399 \to 0.899422073$;
- the transported intervals: $0.001001626 \to 0.001821437$ on $[3.2059998, 3.206]$ and
  $0.001442884 \to 0.002590164$ on $[3.9979773, 4]$;
- the lifts $n_0 = 24, 48, 96$:
  $0.660233530, 0.402885117, 0.063967113 \to 0.660561928, 0.403276789, 0.064969180$;
- the survey seed: $0.543812926 \to 0.543857012$.

The largest $L = 12$ radius fell from $4.07 \times 10^{-6}$ to $1.54 \times 10^{-6}$. The
transported intervals are kept fixed: bisection would now reach $3.2059997$ and $3.9979771$. The
stage-2 case ids `transport-threshold-*` became `transport-*`.

`results/scan_summary.json` and the numerical points of `results/transfer_curve.json` are still
those of the earliest floating-point scan, which used momentum-resolved Lanczos with the removed
local code; `scripts/scan.jl` (now complete spectra from `float_ring_spectrum`) was not rerun over
the full grid. At each revision `scripts/report.jl --certified-only` refreshed the certified curve
points and the plot, and `results/feasibility.json` was regenerated on the current host.

Earlier, the directory-format standardization moved the stage-2 driver to
`scripts/transfer_certify.jl`, and the source-hash policy change made the core hashes provenance
only. Both regenerations left every certified number unchanged.
