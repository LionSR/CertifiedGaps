The dimerized spin-1/2 Heisenberg chain now has independently reviewed finite-spectrum certificates at delta=0.5, 0.3 and 0.25. All changes are confined to `certificates/dimerized_heisenberg/`, the explicitly requested destination.

| delta | n0 (two-spin cells) | beta0 | verified uniform gap | bound / numerical 16-spin gap |
|---|---:|---:|---:|---:|
| 0.5 | 4 | 6 | 0.5326577398571627 | 0.4481193 |
| 0.3 | 4 | 6.5 | 0.17291646480187198 | 0.2131759 |
| 0.25 | 4 | 7.25 | 0.057966591686137485 | 0.0817182 |

These bounds hold for every **even L >= 4 cells**, i.e. **8,12,16,... physical spins**, under the repository's **general internally reflected RP theorem**. The original stricter real identical-factor requirement is false; the exact algebra audit proves that obstruction and does not claim it has been satisfied. The singlet-i basis produces complex conjugate factors. Reference gaps, ratios, unsuccessful scans and correlation-length extrapolations remain numerical.

The proof documents complete dense Sz/momentum/spin-flip sectors, exact Hamiltonian interval construction, rigorous norm bounds for residual and orthogonality errors, polar/Weyl eigenvalue enclosures, outward partition sums, exact rational purity defects and `certified_gap_bound`. Classical binary64 CPU OpenBLAS is an explicit trusted arithmetic assumption. No Lanczos spectra enter certificates.

Validation:

- Original fast suite: **2403/2403**, about 34 seconds; pinned replay on current main: **2403/2403**.
- Independent acceptance: **4012/4012**, 145.5 seconds; regenerated all six full spectra, proved containment inside every stored eigenvalue interval, and checked exact rational residuals in all 63 eight-spin sectors. Floating-point reproducibility differences at delta=0.5 are rigorously covered.
- Current `CertifiedGaps` reconstruction: **1536/1536**, 11.5 seconds, reproducing all three bounds exactly. Validated core tree: `5ebfb34b3f28f37a6fe600857b387dfe887de8b8`.
- No admissible even seed on the tested grid for delta<=0.2. The largest actual ring is 16 spins; no 20- or 40-spin computation is claimed. Odd n0=3 is diagnostic only. Under a 20-spin largest-ring ceiling, the eligible even seed sizes are n0=2 and4. The alternate interpretation n0<=10 would require up to40 spins and remains uncompleted.

Historical source hashes and production results are immutable and pinned to `ef0a3e6cf54ec56dcef5e8a37963d0184f6fd31b`. The current package rename is handled locally: `replay_pinned.jl` extracts that audited core from git into an ignored cache; `validate_migration.jl` separately reconstructs the certificate with current `CertifiedGaps`. No core file or CI is changed.

Reproduce from repository root with Julia 1.12.4 and dependencies:

```sh
OPENBLAS_NUM_THREADS=1 julia --project=julia certificates/dimerized_heisenberg/replay_pinned.jl
OPENBLAS_NUM_THREADS=1 julia --project=julia certificates/dimerized_heisenberg/validate_migration.jl
OPENBLAS_NUM_THREADS=1 julia --project=julia certificates/dimerized_heisenberg/replay_pinned.jl review/independent_review.jl
```

See `REPORT.md`, `ENCLOSURE_PROOF.md`, `review/independent_review.md`, exact-rational JSON results and execution logs. Draft for mathematical and implementation review; do not merge automatically.
