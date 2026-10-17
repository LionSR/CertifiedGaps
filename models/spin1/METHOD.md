# Reducing the computational cost of the spin-1 certificate

The periodic matrix product state (MPS) trial state used in the spin-1 certificate has bond
dimension $D = 26$, parameterized by 23 integers: eight diagonal amplitudes and fifteen couplings
between virtual $\mathrm{SU}(2)$ multiplets. Its virtual Hilbert space decomposes under
$\mathrm{SU}(2)$ as
```math
V_{\mathrm{virt}} = 4V_{1/2} \oplus 3V_{3/2} \oplus V_{5/2}.
```
To evaluate the norm $\langle\psi|\psi\rangle$ and energy $\langle\psi|H_N|\psi\rangle$ on rings
of length $N = 72$ and $N = 120$, one must raise the $D^2 \times D^2 = 676 \times 676$ double-layer
transfer matrix $T$ to high powers. While the upstream reference verifier block-diagonalizes $T$
into eleven $S^z$ charge-difference blocks (of dimension up to 162), we can exploit the full
$\mathrm{SU}(2)$ rotational symmetry of the trial state to reduce $T$ to six total-spin
multiplicity blocks (of dimension at most 56) while preserving every exact integer in the
certificate.

## What the upstream implementation already does

The upstream files linked below are placed in the gitignored `references/` directory by running
`julia tools/fetch_upstream.jl`.

- **Trial contractions:** The pinned
  [upstream trial code](../../references/periodic-gap-source/verification/computations/trials/original_trial.py)
  builds the $26 \times 26$ integer MPS matrices, block-diagonalizes the $676 \times 676$
  transfer matrix by virtual $S^z$ difference, and evaluates matrix powers by repeated squaring.
  Its eleven $S^z$-sector block dimensions are:
  ```
  1, 8, 32, 80, 136, 162, 136, 80, 32, 8, 1.
  ```
- **Thermal trace recursions:** The
  [source-1 thermal runner](../../references/periodic-gap-source/verification/computations/source1/run_thermal.py)
  uses magnetization ($S^z$) blocks, dihedral symmetry representatives, 32-column batches, exact
  fixed-point integer arithmetic, and a global process queue ordered from the largest sectors
  downward. The
  [source-2 thermal runner](../../references/periodic-gap-source/verification/computations/source2/thermal.py)
  packs 192 integer lanes into wide Python integers, including separate coordinates for the
  $2\pi/3$ complex twist.

Although those upstream reductions are already substantial, our initial direct Julia port
(`results/reference/`) left significant performance on the table by relying on generic `BigInt`
matrix multiplication for the trial state, generic 64-bit integer multiplications inside the
thermal inner loops, and an unbalanced thread schedule.

## Verified $\mathrm{SU}(2)$ reduction for MPS trial contractions

In `src/trials_su2.jl`, we construct exact rational virtual raising and lowering operators $J_+$
and $J_-$ directly from the virtual multiplet structure $4V_{1/2} \oplus 3V_{3/2} \oplus V_{5/2}$
and verify all six tensor commutator relations with the physical spin-1 operators by exact
arithmetic. On the $676$-dimensional double layer $V_{\mathrm{virt}} \otimes V_{\mathrm{virt}}^*$,
we form the combined ladder operators
```math
\mathcal{R} = 2J_+ \otimes I - I \otimes J_-, \qquad \mathcal{L} = J_- \otimes I - 2I \otimes J_+.
```
The kernels of $\mathcal{R}$ at each total spin $`J \in \{0, 1, 2, 3, 4, 5\}`$ isolate the
highest-weight multiplicity spaces, which have dimensions
```
26, 56, 48, 24, 7, 1
```
(Notice that $\sum_{J=0}^5 (2J+1) \dim M_J = 1\cdot 26 + 3\cdot 56 + 5\cdot 48 + 7\cdot 24 + 9\cdot 7 + 11\cdot 1 = 676$.)

Rather than relying on numerical eigensolvers, this block decomposition is proved by exact
rational linear algebra. Within **each of the eleven original $S^z$ charge blocks**, applying the
lowering operator $\mathcal{L}$ to the highest-weight basis vectors constructs an explicit square
change-of-basis matrix $W$. The verifier checks in exact arithmetic that $W$ has full rank and
intertwines both the double-layer transfer matrix $T$ and the two-site energy insertion $O$ with
the reduced multiplicity blocks $(T_J, O_J)_{J=0}^5$. These exact intertwining identities
guarantee that for every ring length $N$, the norm and energy contractions equal
```math
V_N = \sum_{J=0}^5 (2J+1)\mathrm{Tr}(T_J^N), \qquad U_N = \sum_{J=0}^5 (2J+1)\mathrm{Tr}(T_J^{N-2}O_J).
```

In addition to computing $V_N$ and $U_N$, the verifier reconstructs the individual trace
contributions of all eleven original $S^z$ charge blocks and checks them against the upstream
length-120 evidence file. Both the length-72 and length-120 numerator and norm integers match the
upstream certificate bit-for-bit.

By working in the $\mathrm{SU}(2)$ multiplicity spaces, the largest matrix dimension raised to the
$`120`$th power shrinks from $162$ to $56$, reducing the $\sum d_i^3$ matrix-multiplication work
proxy from $10{,}373{,}002$ to $317{,}952$ (a factor of $32.6\times$). The remaining exact
integer matrix powers are evaluated using
[Nemo's FLINT-backed integer matrices](https://nemocas.github.io/Nemo.jl/stable/matrix/) (`ZZMatrix`):
rational denominators in $T_J$ and $O_J$ are cleared explicitly, powers are computed over $\mathbb{Z}$,
and every final rescaling division is asserted to have zero remainder. Native integer construction
of the small-entry transfer blocks is also guarded by an explicit overflow bound.

The original `BigInt` charge-block implementation remains available via
`run_trials(...; backend=:reference)` for cross-checking; the default backend is `:su2`.

## Thermal kernel and scheduling optimizations

In the fixed-point thermal trace recursion (`src/thermal.jl`), the off-diagonal bond matrix
elements in the computational basis take only the small integer values $`\{-2, -1, 1, 2\}`$. The
optimized inner-loop kernels expand multiplications by these fixed coefficients into bit shifts,
additions, and subtractions, enabling integer SIMD vectorization without general 64-bit integer
multiplications on every bond hop. Similarly, the source-2 kernel in the Eisenstein basis
$`\{1, \omega\}`$ expands its three fixed coefficient pairs into additions and subtractions, and
applies the arithmetic right-shift floor division to output rows while they remain resident in L1
cache. The fixed-point recurrence formula and a priori rounding-error bounds are completely
unchanged.

Unit tests verify that the optimized kernels produce bit-identical output coordinates to the
reference kernels across edge cases (including signed values near the fixed-point scale and
nonstandard batch widths), and that the full thermal trace sums and rational interval enclosures
match the immutable upstream certificates.

Finally, the multi-threaded task queue orders work globally from the largest sector down across
all chain lengths, boundary twists, and magnetization sectors using Julia's greedy dynamic
scheduler, eliminating the tail-latency imbalance of the initial port.

## Benchmark measurements

The directory `results/reference/` preserves the completed initial Julia reproduction, while
`results/` stores the optimized run. All timings below were measured locally in a single process
to isolate algorithmic and kernel improvements:

- **MPS trial contractions** (`results/trial-benchmark.json`, `scripts/benchmark_trials.jl`):
  Both $N = 72$ and $N = 120$ were contracted by all three methods in a single process after
  JIT warm-up (excluding basis geometry setup equally across methods). All $U_N$, $V_N$, and
  individual charge-block integers matched identically:

  | Exact contraction method | Time for $N = 72$ and $N = 120$ | Speedup vs. baseline |
  | --- | ---: | ---: |
  | Reference Julia `BigInt` ($S^z$ charge blocks, max dim 162) | 33.181536 s | $1.0\times$ |
  | FLINT `ZZMatrix` ($S^z$ charge blocks, max dim 162) | 1.553376 s | $21.4\times$ |
  | FLINT `ZZMatrix` ($\mathrm{SU}(2)$ multiplicity blocks, max dim 56) | 0.092126 s | $360.2\times$ |

  Switching from generic `BigInt` arrays to FLINT `ZZMatrix` yields a $21.4\times$ speedup, and
  the $\mathrm{SU}(2)$ multiplicity reduction contributes an additional $16.9\times$ speedup on
  top of FLINT.
- **Thermal trace kernels** (`results/kernel-benchmark.json`, `scripts/benchmark_kernels.jl`):
  Warmed single-sector contractions in one process achieved a $2.18\times$ speedup for source 1
  (length 11, $\pi$ twist, magnetization $M = 1$) and a $2.58\times$ speedup for source 2
  (length 10, $2\pi/3$ twist, magnetization $M = 1$), with identical integer outputs.
- **Column batch width** (`results/batch-benchmark.json`, `scripts/benchmark_batches.jl`):
  Evaluating the same 64 columns of the largest real-twist sector at batch widths $4, 8, 16,$ and
  $32$ produced identical integers; smaller batch widths were slower, confirming $32$ as the
  default.
- **End-to-end component wall times:** See [`README.md`](README.md#reproduce) for the full
  4-thread reproduction timings recorded in `results/`.
