# Reducing the cost of the spin-1 certificate

The physical trial has bond dimension 26, but only 23 integer parameters: eight diagonal amplitudes and fifteen couplings between virtual multiplets. Its virtual representation is

\[
4V_{1/2}\oplus3V_{3/2}\oplus V_{5/2}.
\]

The source verifier contracts its 676-dimensional double layer in eleven charge-difference blocks. We can use the full spin symmetry to reduce that calculation further while preserving every integer in the certificate.

## What the upstream implementation already does

The upstream files linked here live in the gitignored `references/` directory; `julia tools/fetch_upstream.jl` fetches them.

The pinned [upstream trial code](../../references/periodic-gap-source/verification/computations/trials/original_trial.py) builds the integer matrices, separates charge differences, and evaluates powers by repeated squaring. Its block sizes are

```
1, 8, 32, 80, 136, 162, 136, 80, 32, 8, 1.
```

The [source-1 thermal runner](../../references/periodic-gap-source/verification/computations/source1/run_thermal.py) uses magnetization blocks, dihedral representatives, 32-column batches, exact fixed-point arithmetic and a global process queue starting with the largest chains. The [source-2 calculation](../../references/periodic-gap-source/verification/computations/source2/thermal.py) packs 192 integer lanes into large Python integers, including separate coordinates for the complex twist.

These are substantial existing reductions. The first Julia reproduction preserved the arithmetic but left performance on the table: generic `BigInt` matrix multiplication, generic integer multiplications in the thermal inner loops, and an initially unbalanced sector schedule.

## Verified SU(2) reduction

`src/trials_su2.jl` constructs rational virtual raising and lowering matrices directly from the multiplet labels. All six tensor commutator relations are verified exactly. On the double layer it uses

\[
\mathcal R=2J_+\otimes I-I\otimes J_-,\qquad
\mathcal L=J_-\otimes I-2I\otimes J_+.
\]

The highest-weight kernels of `R` give multiplicity spaces for total spin 0 through 5, of dimensions

```
26, 56, 48, 24, 7, 1.
```

This is checked by exact rational algebra, not inferred from numerically clustered eigenvalues. For **each of the eleven original charge blocks**, lowering the highest-weight vectors constructs a square matrix `W`. The code verifies that `W` has full rank and intertwines both the transfer matrix `T` and energy insertion `O` with the reduced blocks. These identities establish an exact change of basis for both matrices and therefore preserve the contractions at every length.

The contractions become

\[
V_N=\sum_{J=0}^5(2J+1)\operatorname{Tr}(T_J^N),\qquad
U_N=\sum_{J=0}^5(2J+1)\operatorname{Tr}(T_J^{N-2}O_J).
\]

The program also reconstructs the contribution of every original charge block and compares the full length-120 block integers with upstream evidence. Both full length-72 and length-120 numerator/norm integers match exactly.

The largest powered matrix shrinks from 162 to 56. The sum-of-cubes matrix-size proxy falls from 10,373,002 to 317,952, a factor of 32.6; this is a work proxy, not a measured runtime ratio. [Nemo’s FLINT-backed integer matrices](https://nemocas.github.io/Nemo.jl/stable/matrix/) perform the remaining exact matrix powers. Rational block denominators are cleared explicitly, powers are taken over the integers, and every final division is checked for zero remainder. Native integer construction of the small-entry transfer blocks also has an explicit overflow bound.

The original `BigInt` charge-block implementation remains available as `run_trials(...; backend=:reference)` for comparison. The default is `:su2`.

## Thermal changes

The optimized kernels expose the exact off-diagonal coefficients `-2`, `2`, `-1`, and `1` as shifts, additions and subtractions. This permits integer SIMD without general 64-bit multiplications for every hop. The Eisenstein-coordinate kernel similarly expands its three fixed coefficient pairs into additions and subtractions. Output rows are floor-divided while still in cache. The error bounds and integer recurrence are unchanged.

Tests compare every output coordinate with the reference kernels, including signed values near the fixed-point scale and a nonstandard batch width. Full thermal sums and rational enclosures are then compared with the immutable source certificates.

The task queue now follows the upstream largest-first approach across all lengths, twists and magnetization sectors, with Julia's greedy scheduling. This corrects the unbalanced scheduling of the first port. Sector completion is reported during long runs.

## Measurements

`results/reference/` preserves the completed first Julia reproduction and `results/` the optimized run. Every measurement below is local and scoped as stated; none is a comparison with the upstream authors' hardware.

- **Trials** (`results/trial-benchmark.json`, `scripts/benchmark_trials.jl`). Both N=72 and N=120 were contracted by all three methods in one process, after warming each; geometry construction was excluded equally, and the complete U/V integers and every charge-block contribution agreed.

  | Exact contraction method | Time for both lengths |
  | --- | ---: |
  | Original Julia `BigInt` charge blocks | 33.181536 s |
  | FLINT with the same charge blocks | 1.553376 s |
  | FLINT with SU(2) multiplicity blocks | 0.092126 s |

  The warmed contraction speedup over the first Julia implementation is 360.2x, of which the SU(2) reduction contributes 16.9x beyond the change of arithmetic backend.
- **Thermal kernels** (`results/kernel-benchmark.json`, `scripts/benchmark_kernels.jl`). Warmed reference and optimized sector contractions in one process gained about 2.18x for source 1 (length 11, twist pi, magnetization 1) and 2.58x for source 2 (length 10, complex twist, magnetization 1). The complete integers agree.
- **Batch width** (`results/batch-benchmark.json`, `scripts/benchmark_batches.jl`). The same 64 columns of the largest real-twist sector at batch sizes 4, 8, 16 and 32 give identical integers; smaller batches were slower, so the default stays at 32. This is a partial-column experiment, not a full thermal certificate.
- **Full components.** The `wall_seconds` recorded in the component outputs include cold compilation and any competing work; see the README for the stored run.

The physical trial remains bond dimension 26; reducing it would need a new variational search and a new check of the very tight energy threshold. The same representation-theoretic reduction is a natural starting point for an SU(2)-structured spin-2 trial, with its own tensor relations and finite checks.
