# Spin-1 periodic certificate reproduction in Julia

This is the reference calculation for extending the thermal-MPO gap method to spin two. It recomputes the finite certificates of **The periodic spin-one Haldane gap** (OpenAI, September 24, 2026) from their integer recipes. ITensors supplies independent physical-MPO, MPS, DMRG, and thermal-MPO controls.

The certification path uses exact integers and rational error bounds. It does not infer a certificate from a converged DMRG energy or a truncated floating-point MPO.

## Reproduce

From the repository root, with Julia 1.12:

```sh
julia --project=models/spin1 -e 'using Pkg; Pkg.instantiate()'
julia --startup-file=no --project=models/spin1 models/spin1/test/runtests.jl
julia --startup-file=no --threads=4 --project=models/spin1 models/spin1/scripts/reproduce.jl all /tmp/spin1-fresh-run
```

`spin1` is the package `Spin1Certificates`; every script, test and experiment loads it as a package (`using Spin1Certificates` or `import Spin1Certificates as SC`) instead of including `src/`. It depends on the shared core [`CertifiedGaps`](../julia) as a path dependency, so instantiate from a checkout that contains `julia/`. The core supplies the exact-arithmetic and result helpers (`RationalBig`, `frac`, `parse_rational`, `parse_exact_json`, `write_result`, `source_hashes`), the rational doubling chain with the bootstrap tail (`Recurrence.rational_gap_chain`, `squaring_defect`, `bootstrap_coefficient`) and the ITensors thermal purification used by the numerical control (`Purification.purified_log_partition`). The ITensors controls live in the submodule `Spin1Certificates.ITensorsChecks`. The implementation fingerprint in `verification.json` covers spin1's own environment and sources, not the core's files.

Choose an output directory for each fresh run. The default is `models/spin1/results`. The complete runner recomputes every contraction, performs the ITensors controls, and writes `verification.json` only as passed after every component succeeds. The manifest pins ITensors 0.9.32, ITensorMPS 0.4.2, JSON3 1.14.3, Nemo 0.57.0 and the other dependencies. No Python numerical packages are needed to run the Julia port.

To check the included outputs without repeating the expensive contractions:

```sh
julia --startup-file=no --project=models/spin1 models/spin1/scripts/reproduce.jl audit models/spin1/results
julia --startup-file=no --project=models/spin1 models/spin1/test/audit_failures.jl models/spin1/results
```

The audit explicitly reports `execution_mode: stored-output-audit`. It checks all full contraction integers against immutable upstream evidence, rebuilds their rational intervals, checks sector sums, and recomputes the inexpensive arithmetic and scalar assertions. It does **not** claim a fresh thermal or long-chain trial contraction. A failed audit replaces any earlier success record with `status: failed`.

The components can also be run separately:

```sh
julia --startup-file=no --threads=4 --project=models/spin1 models/spin1/scripts/thermal.jl 1
julia --startup-file=no --threads=4 --project=models/spin1 models/spin1/scripts/thermal.jl 2
julia --startup-file=no --project=models/spin1 models/spin1/scripts/trials.jl
julia --startup-file=no --project=models/spin1 models/spin1/scripts/scalars.jl
julia --startup-file=no --project=models/spin1 models/spin1/scripts/arithmetic.jl
julia --startup-file=no --project=models/spin1 models/spin1/scripts/itensors_crosscheck.jl
```

Both thermal scripts accept an optional output filename after the source number; the other component scripts accept an output filename as their first argument. The largest computation is source 1 at length 12. Magnetization sectors across all lengths and twists run on the available Julia threads in a global, largest-first greedy queue. Work within each sector is batched in 32 integer columns.

## What is checked

We use the periodic Hamiltonian

\[
H_n=\sum_{j=1}^{n}\mathbf S_j\cdot\mathbf S_{j+1},\qquad
\mathbf S^2=2I,\quad S_{n+1}=S_1,
\]

with bond coupling one and shifted partition functions
\(Z_n(\beta,g)=\operatorname{Tr}\exp[-\beta(H_n^g+naI)]\), where \(a=700741/500000\). The twists are identity, a rotation by \(\pi\), and a rotation by \(2\pi/3\) at the closing bond. The inverse-temperature shift cancels from both purity ratios.

| Component | Fresh computation | Comparison |
| --- | --- | --- |
| `src/thermal.jl`, source 1 | 18 traces: lengths 4–12, identity and pi twists, beta 49/4 | Every full integer total and rational interval endpoint |
| `src/thermal.jl`, source 2 | 9 traces: lengths 6, 8, 10, all three twists, beta 21/2 | Every full integer total, error bound, and rational interval endpoint |
| `src/trials.jl` | Rebuild all 2,028 entries of three 26 by 26 integer matrices; contract periodic MPS at lengths 72 and 120 | Full energy numerators, norms, positive variational margins; all length-120 transfer blocks |
| `src/arithmetic.jl` | Independent column-action Hamiltonians, all magnetization blocks, twisted symmetries, 27 weighted triad rows | Exact equality and integer/rational error bounds |
| `src/scalars.jl` | Moment filters, capped masses, both purity initializations, six rounded updates, bootstrap constants | Rational inequalities and every published rounded pair |
| `src/itensors_checks.jl` (`ITensorsChecks`) | Four-site MPOs for all twists, integer trial imported as MPS, QN DMRG, thermal MPO | Independent physical matrices, exact four-site energy, and Trotter step refinement |

The integer thermal recurrence uses arithmetic right shifts, which implement floor division also for negative entries. Explicit bounds rule out overflow in `Int64` partial sums and `Int128` norm chunks; final totals and trial contractions use `BigInt`. Source 2 is an independent unpacked implementation in the basis \(1,\omega\), with \(\omega^2+\omega+1=0\), rather than a translation of the upstream packed-lane bit tricks. JSON integer tokens are preserved before parsing, so large certificate integers never pass through `Float64`.

The default trial verifier uses exact SU(2) multiplicity blocks of sizes 26, 56, 48, 24, 7 and 1, with a rational change of basis verified against every original charge block. Nemo/FLINT performs the integer powers. See [`OPTIMIZATION.md`](OPTIMIZATION.md) for the derivation, upstream comparison and benchmarks.

The exact trials establish

\[
\langle H_{72}\rangle/72<-1.401482,
\qquad \langle H_{120}\rangle/120<-1.401482.
\]

For the two-purity definitions
\(S(n,\beta)=Z_{2n}(\beta)/Z_n(\beta)^2\) and
\(T(2n,\beta)=Z_{2n}(2\beta)/Z_{2n}(\beta)^2\), the reproduced finite arithmetic supplies:

| Initialization | Spatial purity | Physical purity | Parameters for the paper's bootstrap theorem |
| --- | --- | --- | --- |
| Source 1 | \(S(2304,784)\ge 1-84513/10^7>0.99\) | \(T(4608,784)\ge 1-27493/(5\cdot10^6)>0.99\) | \(u_0=1/100,\ C=5\) |
| Source 2 | \(S(60,105/4)>7/8\) | \(T(120,105/4)>7/8\) | \(u_0=1/8,\ C=79/10\) |

Combined with the analytic transfer and bootstrap arguments in the paper, these are the finite inputs for gap bounds \(\log(20)/784\) for even lengths at least 2304 and \((4/105)\log(80/79)\) for even lengths at least 60. The Julia program verifies finite arithmetic; it does not formally prove the transfer construction or the all-scale induction. The boundary-field companion and a spin-2 certificate are outside this reproduction.

## Stronger bounds from the same finite inputs

The new [temperature-lifting calculation](experiments/temperature_lift.jl) establishes **gap > 0.072901771 for every even length at least 36**, combining its exact finite checks with the derivation in the note and the original analytic transfer/bootstrap theorems. Simply doubling inverse temperature with a length-30 moment bound proves **gap > 7/100 = 0.07**. It combines the eighteen source-1 thermal values with the existing length-120 trial, raises inverse temperature before extracting purity bounds, and retains the threefold sector multiplicity when bounding residual moments. The earlier length-72 seed covers lengths 36–72. It needs no new thermal or trial contractions, and no larger MPS.

```sh
julia --startup-file=no --project=models/spin1 models/spin1/experiments/temperature_lift.jl
julia --startup-file=no --project=models/spin1 models/spin1/experiments/test_temperature_lift.jl
```

The [rational certificate](experiments/temperature-lift.json) records all fractional-power enclosures, dominant-weight refinements, bootstrap steps, and the bridge covering the early lengths. This is a new analytic derivation with exact finite arithmetic, not an independently reviewed or formally verified theorem.

For reproducible numerical discovery, `julia --startup-file=no --project=models/spin1 models/spin1/experiments/search_temperature_lift.jl` scans the moment length and temperature multiplier. Its output is explicitly noncertifying; the selected parameters must pass the rational verifier above. The new tests passed 285 checks, including forged-certificate rejection, and all 549 earlier scalar checks still pass.

The separate [scalar optimization](experiments/optimize_bounds.jl) tightens the existing polynomial cutoffs, retains unequal purity defects during the coupled updates, and uses the interpolation inequality to cover the early lengths before invoking the uniform bootstrap theorem. Its exact rational checks and the paper's analytic arguments give **gap > 0.016978269 for every even length at least 36**. A compact version with six-decimal cutoffs gives **gap > 1/60** in the same range. The original reproduction outputs remain the baseline.

Using only the less expensive source-2 inputs gives **gap > 0.011977901 for every even length at least 60**. This uses the same nine short-chain traces and length-120 trial as the original source-2 argument. No new thermal or trial contractions are needed for either scalar improvement.

```sh
julia --startup-file=no --project=models/spin1 models/spin1/experiments/optimize_bounds.jl
julia --startup-file=no --project=models/spin1 models/spin1/experiments/test_bounds.jl
```

See the [derivation](experiments/optimized-spin1-bounds.tex) and [complete rational certificate](experiments/bounds.json). The experiment verifies the hashes of the completed baseline reproduction, so run the baseline reproduction first if those outputs are absent. The final inequalities use rational logarithm enclosures; floating-point values are display aids only. This improves the constants and length coverage, without claiming an optimal bound or formalizing the paper's analytic theorems.

## Provenance and outputs

The unchanged upstream source and evidence are read from `../../references/periodic-gap-source`, pinned to commit `adc7f1241b42e322a6451854ab7e4b4c146bf78a`; they are not stored in this repository. `julia tools/fetch_upstream.jl` (from the repository root) downloads the files below at that commit and checks the pinned manifest hash and every file hash. The original download checked each Git blob hash and recorded SHA-256 hashes in `download-provenance.json`. Every Julia component pins the upstream source-manifest hash and checks 43 source/evidence files before running.

See [`NOTICE`](NOTICE) and the preserved upstream Apache-2.0 license. These are modified Julia implementations, not an upstream Julia release.

[`results/verification.json`](results/verification.json) records component hashes and the current verifier/environment fingerprints. The initial implementation outputs are retained in `results/reference/`; the primary outputs are from the complete optimized rerun. The component JSON files retain full exact integers, rational bounds, run timestamps and timings. Rationals use Julia's `numerator//denominator` string notation. Numerical ITensors outputs are marked `certified: false`. The reproduction report is [`../working/reproduction/reproduction-report.md`](../working/reproduction/reproduction-report.md).
