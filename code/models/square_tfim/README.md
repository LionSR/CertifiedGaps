# Square-lattice transverse-field Ising model

An analytic theorem gives a unique ground state and a uniform spectral gap for the
ferromagnetic transverse-field Ising model on every finite open or periodic rectangle when
$g>4J$. The proof sums Ising paths in imaginary time and uses FKG to control the full
spectral gap. Independent spacetime bootstrap certificates at $g/J=18,20,24$ are retained.

| parameter $g/J$ | seed $(n_0=m_0,\beta_0 J)$ | certified gap $\ge$ | lengths covered |
|---:|---|---:|---|
| $>4$ | analytic theorem | $2\sqrt{g(g-4J)}$ | all finite rectangles |
| 4.1 | analytic theorem | $1.28062J$ | all finite rectangles |
| 4.25 | analytic theorem | $2.06155J$ | all finite rectangles |
| 4.5 | analytic theorem | $3J$ | all finite rectangles |
| 5 | analytic theorem | $4.47213J$ | all finite rectangles |
| 6 | analytic theorem | $6.92820J$ | all finite rectangles |
| 8 | analytic theorem | $11.31370J$ | all finite rectangles |

## Model

With Pauli matrices and periodic indices, the Hamiltonian is

```math
H_{n,m}=-J\sum_{x,y}(Z_{x,y}Z_{x+1,y}+Z_{x,y}Z_{x,y+1})-g\sum_{x,y}X_{x,y}.
```

The computation sets $J=1$; energies and inverse temperatures rescale as in the table.
Lengths count individual lattice spacings, and an $L$ by $L$ torus has $L^2$ spins. At period
two each bond in that direction is counted twice, as required by the periodic sum and the
transfer trace. The bounds also cover $L=2$ with this convention; for $L\ge3$ it is the
ordinary square lattice with each undirected bond counted once. The analytic result covers
open boundaries by omitting the boundary-crossing bonds.

## Reflection positivity

In the physical basis, a crossing bond has $-h=Z\otimes Z$. After a Hadamard rotation on
every spin, a complete row has inter-row bond $-\sum_y X_y\otimes X_y$ and real symmetric
on-site term $-\sum_y X_yX_{y+1}-g\sum_y Z_y$. The integer factors $C_y=X_y$ are checked by
`SquareIsing.row_rp` and the core `verify_rp`.

In both spatial directions, finite Trotter transfers have the form $D^{1/2}KD^{1/2}$ with
real symmetric positive K. Their scalar squaring and moment interpolation inequalities
pass to the exact quantum partition function. See [METHOD.md](METHOD.md) for the proof;
no numerical Trotter or tensor truncation is used in the certificate.

## Seeds and certified bounds

The general graph theorem assumes nonnegative pair couplings, fields $g_i\ge h>0$, and an
incident coupling bound $\sum_jJ_{ij}\le s<h$. It proves a unique ground state and
$\Delta\ge2\sqrt{h(h-s)}$, uniformly in graph size. Taking $h=g$ and $s=4J$ gives the
headline table. The complete proof is in [METHOD.md](METHOD.md); exact rational evaluations
are in `results/high_field_bounds.json`. No finite purity seed enters this argument.

The independent spacetime bootstrap records are:

| $g/J$ | $(n_0=m_0,\beta_0J)$ | bootstrap gap $\ge$ |
|---:|---|---:|
| 18 | $(2,53/240)$ | $0.74378J$ |
| 20 | $(2,157/800)$ | $2.51376J$ |
| 24 | $(2,79/480)$ | $5.52173J$ |

These points use the staggered defects of $P_x(2,2,\beta_0)$,
$P_y(4,2,\beta_0)$ and $P_\tau(4,4,\beta_0)$. The four log partition enclosures come from
complete verified spectra of the 2 by 2, 4 by 2, and 4 by 4 tori. The largest symmetry block
has dimension 2392. The exact defects and every rational recurrence step are in
`results/certificate_g18.json.gz`, `certificate_g20.json.gz`, and `certificate_g24.json.gz`.

The refined seed bounds are exact terminating decimals:

| $g/J$ | $\varepsilon_x$ | $\varepsilon_y$ | $\varepsilon_\tau$ |
|---:|---:|---:|---:|
| 18 | 0.02456152 | 0.04825847 | 0.01222144 |
| 20 | 0.01967145 | 0.03878343 | 0.01325307 |
| 24 | 0.01375692 | 0.02723951 | 0.01236053 |

The code uses the three-component update from the
[2+1 dimensional note](../../../paper/thermal-mpo-gap-2plus1d.tex), then its quadratic tail
after twelve steps, when all defects are at most $1/1000$. A rational grid of denominator
$10^{2000}$ resolves the defects at these scales. Exact lower bounds for the finitely many
earlier interpolation windows ensure that the conclusion covers every even square side,
not only a sequence of doubled tori. The one-dimensional `certified_gap_bound` is not used.

The `refinement` field contains the new seeds and bounds; the original `seed` and `chain`
fields retain their previous exact values. The checker recomputes both. Those records also
contain verified finite-volume ground energies and gaps.

At smaller fields, complete **floating-point** diagonalization gives the following
sixteen-spin values. These are numerical diagnostics, separate from the analytic bound.

| $g/J$ | $E_0(4,4)/(16J)$, approximately | $\Delta_{4,4}/J$, approximately |
|---:|---:|---:|
| 4 | −4.13897624 | 3.24828451 |
| 5 | −5.10617353 | 5.44369104 |
| 6 | −6.08666146 | 7.56254563 |
| 8 | −8.06381620 | 11.69479751 |

No thermodynamic extrapolation is used. At these four fields none of the tested
$2\times2$ seeds closes on the grid $\beta g=0.5,0.525,\ldots,16$. At $g/J=4$ the analytic
bound also vanishes; the finite-volume gap does not prove a uniform positive gap there.

The floating-point search at $g/J=17$ found no closing seed on the tested grid
$\beta g=2,2+1/8,\ldots,12$. This does not exclude another temperature, a larger seed, or
a stronger recurrence, and says nothing about the physical critical field.

## Reproduce

From the repository root, using Julia 1.12 and the core environment:

```sh
export OPENBLAS_NUM_THREADS=1
julia --project=code/julia code/models/square_tfim/certify.jl
julia --project=code/julia code/models/square_tfim/check.jl
julia --project=code/julia code/models/square_tfim/check.jl --full
julia --project=code/julia code/models/square_tfim/test/runtests.jl
julia --project=code/julia code/models/square_tfim/scripts/high_field.jl --check
```

The default drivers process the three bootstrap fields and the analytic arithmetic record.
The last command checks only the analytic record, without diagonalization. Omit `--check`
to regenerate it. A bootstrap field can be supplied, for example
`check.jl 18 --full`. `certify.jl --reuse` reuses the stored verified spectra after checking
their hashes and row models; the default regenerates them. On the development Apple Silicon
machine with one BLAS thread, the three refined certificates regenerated from stored spectra
in 26.2 seconds, and the 58 model tests took 133.3 seconds in their test sets. The original
compressed fast checker took 38 seconds; full spectral generation took 66--70 seconds per
field. Full spectral replay of all three refined certificates took 217.1 seconds and passed.
Measurements exclude startup and initial package precompilation. Certificates and
spectral records use lossless gzip JSON.

The additional 25 high-field tests take 1.5 seconds after loading. They include independent
diagonalizations of eighteen weighted paths, stars, and complete graphs with nonuniform
fields. The analytic-only checker takes 2.35 seconds including Julia startup.

For numerical exploration only:

```sh
julia --project=code/julia code/models/square_tfim/scripts/scan.jl 17
julia --project=code/julia code/models/square_tfim/scripts/refine_scan.jl
julia --project=code/julia code/models/square_tfim/scripts/moderate_fields.jl
```

The first command takes about 35 seconds per field. The second reuses stored spectra to
search temperature with the longer recurrence. The third computes fields 4, 5, 6, and 8,
with about 25 seconds per sixteen-spin spectrum. All three write floating-point diagnostics.

The general Lean theorem can be checked with
`cd code/lean && lake build CertifiedGaps.CubeCriterion`.

## What is rigorous

The graph high-field theorem is a mathematical proof using the classical FKG inequality,
the positive Ising path expansion, and a finite-dimensional Trotter limit. Its numerical
evaluations use only exact integer and rational arithmetic, including an integer square
root. `check_high_field` checks those evaluations and source hashes; it does not check the
mathematical proof. This graph theorem has not been formalized in Lean.

For the bootstrap certificates, the spatial transfer argument, exact RP identity,
full spectral enclosures, interval
partition sums, and exact rational three-direction recurrence establish the bounds together.
The spectral arithmetic uses the shared `VerifiedSpectra` implementation. Its trusted base
is IEEE binary64 round-to-nearest with gradual underflow, the standard a priori dot-product
rounding bound for any summation order and FMA, BLAS/LAPACK finite outputs, and
`IntervalArithmetic`; residual and orthogonality bounds are followed by Weyl and
Ostrowski/polar-factor estimates. All logarithms in the final gap bound use the exact
rational series in `ExactArithmetic.log_enclosure`.

The fast checker replays stored spectral enclosures and checks their hashes. `--full` rebuilds
every symmetry block and independently verifies that the fresh enclosures imply the stored
seed bounds. The tests compare small tori with a separately constructed Hamiltonian, check
the full sixteen-spin spectral moments, and reject altered seeds, gaps, and spectrum hashes.
The seed scans are floating-point diagnostics only. The Lean theorem
`CertifiedGaps.Spacetime.gapGT_of_cube_seeds` proves the general small-defect criterion from
scalar moment and spatial interpolation hypotheses, without `sorry` or additional axioms.
The Ising transfer construction, the finite asymmetric recurrence preceding the small-defect
regime, and the numerical seeds remain outside that formalization. The model certificate
therefore combines a mathematical proof with verified Julia arithmetic. Its scope is even
square tori with periodic boundaries. The graph theorem has the broader scope stated above.

## Layout

- `src/SquareIsing.jl`: exact row models, transverse symmetries, and RP factors.
- `src/CubeRecurrence.jl`: exact three-component recurrence and interpolation-window bounds.
- `src/Certificate.jl`, `certify.jl`, `check.jl`: certificate generation and replay.
- `src/HighField.jl`, `scripts/high_field.jl`: exact arithmetic for the graph theorem.
- `scripts/build_spectra.jl`, `scripts/scan.jl`, `scripts/refine_scan.jl`: spectral generation
  and numerical seed searches.
- `scripts/moderate_fields.jl`: finite-volume spectra and seed obstructions at fields 4--8.
- `test/runtests.jl`: independent geometry, arithmetic, spectral-moment, and rejection tests.
- `METHOD.md`, `IDEAS.md`: proof and route analysis.
- `results/`: relative symlink to `../../../data/models/square_tfim/`, containing certificates,
  spectrum manifests, compressed verified spectra, and numerical scans.

Model-local `src/`, `scripts/`, `certify.jl`, and `check.jl` hashes are enforced. Core hashes
are provenance only; the checker reruns core arithmetic. Spectrum files are hash-pinned by
the certificates. All can be regenerated in a few minutes; no immutable source exceptions
are needed.
