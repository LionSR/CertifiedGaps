# Audit outcome

**PASS under the repository two-purity theorem and the stated IEEE/BLAS arithmetic
model.** The independent model-specific proof/code audit found no defect. Fresh
full dense spectral replay independently repeated the computation from rational
model inputs and reproduced both original certificates exactly. This is a new
execution of the reviewed implementation, not a second independently written
large-matrix eigensolver or a formal verification of the runtime and theorem.

The audited current main revision is
`22b2759e0790197b83e2571b465aa8dc69ebc0be` (`CertifiedGaps`), verified against
upstream main before and after the replay. Original computation records and
`provenance/legacy/` remain unchanged. New evidence is confined to `audit/`.

| D | Exact seed (n0,beta) | Fresh raw interval-certified lower bound | Published bound, rounded down |
|---|---|---:|---:|
| 4 | (6,9//4) | 0.6946941929576096 | 0.694694192 |
| 1603//500 = 3.206 | (6,381//125) | 0.0022134568960786357 | 0.002213456 |

Both conclusions cover **every even L>=6 one-site cells = sites**, with a unique
ground state, under the repository theorem. D=3.206 is the smallest successful
sampled value; no global threshold over real D is established.

## What was actually recomputed

The historical/pinned replay and the original fast/current checker trusted
stored spectral radii. They reconstructed partition intervals and recurrence
bounds, but did not recompute dense eigenvectors, residuals or orthogonality
errors. They were not full spectral replays.

The new `replay_worker.jl` reads no historical spectral data. For each D it
rebuilds all L=6 and L=12 symmetry blocks from the exact model, densely
recomputes every eigensystem, and reevaluates interval construction errors and
rigorous residual/orthogonality norm bounds. The two charge groups are disjoint.
The assembler independently verifies 29 L=6 blocks (729 levels) and 92 L=12
blocks (531441 levels after multiplicities), source hashes, worker manifests,
Julia/BLAS configuration and unchanged core files. It calculates the fresh
rational defects and recurrence bound before reading historical data for
comparison. No grid search was repeated.

Every sorted eigenvalue center and every radius is bitwise equal to its original
counterpart, at both lengths and both D values. Both exact rational seeds and
raw recurrence bounds also reproduce exactly. The separate certificate JSONs
store the new source/core/input hashes and package metadata:

- [D=4 fresh certificate](results/certificate_D4_1.json)
- [D=3.206 fresh certificate](results/certificate_D1603_500.json)

Per-block spectra and fresh construction/residual/orthogonality evidence are in
`results/D4_1/` and `results/D1603_500/`. Worker source hashes identify the exact
code that executed. The first metadata-type-error launch is retained separately
under `initial-driver-failure/`; no L=12 data from it was accepted.

## Actual validation and runtimes

| Check or worker | Result | Runtime |
|---|---|---:|
| Independent exact Cartesian-to-normalized-Sz bridge | 129/129 | 4.5 s test execution |
| Independent 256-bit residual/orthogonality diagnostics | 90/90 | 5.4 s test execution |
| Original checker on current core | 134/134 | 22.6 s test execution |
| Fresh evidence/rational certificate checker | 240/240 | 17.896 s wall, 10.7 s tests |
| Fresh certificate assembly | Both reproduced exactly | 18.956 s wall |
| D=4 heavy worker | Complete | 1489.782 s = 24.83 min |
| D=4 mixed worker, including L=6 | Complete | 1265.312 s = 21.09 min |
| D=3.206 heavy worker | Complete | 1300.990 s = 21.68 min |
| D=3.206 mixed worker, including L=6 | Complete | 1429.762 s = 23.83 min |

Worker runtimes are measured inside the driver, excluding Julia startup. All
four processes used one BLAS thread; each remained below 25 minutes. Their
logs, manifests and exact timings are retained. The four independent processes
shared the four-CPU host. The fast audit checker verifies this new stored
evidence; it does not itself repeat dense ED. [README.md](README.md) gives the
commands for both levels of replay.

## Scope and limits

[INDEPENDENT_REVIEW.md](INDEPENDENT_REVIEW.md) audits exact RP and normalization,
sector completeness, all three spectral error sources, the polar/Weyl theorem,
the rational purity directions, even-length interpolation and uniqueness.
The 256-bit checks are diagnostics supporting the audited norm inequalities,
not the source of certification.

The trusted base remains Julia, IntervalArithmetic outward rounding, the stated
standard IEEE/BLAS operation model, and the repository moment/recurrence
theorem. Neither the BLAS binary nor the entire theorem has been formally
verified by this audit. The fresh large computations reuse the audited model
and enclosure implementation; the independent reviewer supplies the separate
proof/code review and exact normalization checks. Failed search points, ED
reference gaps/ratios, extrapolations, and statements about all smaller real D
remain numerical or unproved. No model-specific correction or unresolved
replay failure remains.
