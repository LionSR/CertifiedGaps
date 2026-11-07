# Square-lattice transverse-field Ising model

The spacetime bootstrap certifies a unique ground state and a uniform spectral gap for the
ferromagnetic transverse-field Ising model on every even square torus of side at least four.
Both spatial periods grow: these are two-dimensional bounds, not fixed-width ladder bounds.
The proof uses exact spatial transfer inequalities and verified spectra of at most sixteen
spins. The two-dimensional theorem is proved mathematically, but is not yet formalized in Lean.

| parameter $g/J$ | seed $(n_0=m_0,\beta_0 J)$ | certified gap $\ge$ | lengths covered |
|---:|---|---:|---|
| 18 | $(2,2/9)$ | $0.59437J$ | even square side $L\ge4$ |
| 20 | $(2,1/5)$ | $2.17363J$ | even square side $L\ge4$ |
| 24 | $(2,1/6)$ | $4.80557J$ | even square side $L\ge4$ |

## Model

With Pauli matrices and periodic indices, the Hamiltonian is

```math
H_{n,m}=-J\sum_{x,y}(Z_{x,y}Z_{x+1,y}+Z_{x,y}Z_{x,y+1})-g\sum_{x,y}X_{x,y}.
```

The computation sets $J=1$; energies and inverse temperatures rescale as in the table.
Lengths count individual lattice spacings, and an $L$ by $L$ torus has $L^2$ spins. At period
two each bond in that direction is counted twice, as required by the periodic sum and the
transfer trace. The certificate also covers $L=2$ with this convention; for $L\ge4$ it is the
ordinary square lattice with each undirected bond counted once.

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

All points use the staggered defects of $P_x(2,2,\beta_0)$,
$P_y(4,2,\beta_0)$ and $P_\tau(4,4,\beta_0)$. The four log partition enclosures come from
complete verified spectra of the 2 by 2, 4 by 2, and 4 by 4 tori. The largest symmetry block
has dimension 2392. The exact defects and every rational recurrence step are in
`results/certificate_g18.json.gz`, `certificate_g20.json.gz`, and `certificate_g24.json.gz`.

The certified defect bounds are exact terminating decimals:

| $g/J$ | $\varepsilon_x$ | $\varepsilon_y$ | $\varepsilon_\tau$ |
|---:|---:|---:|---:|
| 18 | 0.02470865 | 0.04854442 | 0.01164265 |
| 20 | 0.02002999 | 0.03948414 | 0.01145037 |
| 24 | 0.01392379 | 0.02756791 | 0.01120587 |

The code uses the three-component update from the
[2+1 dimensional note](../../../paper/thermal-mpo-gap-2plus1d.tex), then its quadratic tail
once all defects are at most $1/1000$. Exact rational lower bounds for the finitely many
earlier interpolation windows ensure that the conclusion covers every even square side,
not only a sequence of doubled tori. The one-dimensional `certified_gap_bound` is not used.

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
```

The default drivers process all three fields. A field can be supplied, for example
`check.jl 18 --full`. `certify.jl --reuse` reuses the stored verified spectra after checking
their hashes and row models; the default regenerates them. Measured on the development
Apple Silicon machine with one BLAS thread: spectral generation took 66--70 seconds per
field; `certify.jl --reuse` took 26 seconds for all three fields; the compressed fast checker
took 38 seconds; and full spectral replay took 235 seconds for all three fields. The 42 model
tests took 101--121 seconds in their test sets. These measurements exclude process startup and initial
package precompilation. Certificates and spectral records are stored as lossless gzip JSON;
compression was checked to preserve every seed, recurrence step, and spectral record.

For numerical exploration only:

```sh
julia --project=code/julia code/models/square_tfim/scripts/scan.jl 17
```

This takes about 35 seconds per field and writes a clearly labeled diagnostic, not a certificate.

## What is rigorous

The spatial transfer argument, exact RP identity, full spectral enclosures, interval
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
The seed scans are floating-point diagnostics only. No infinite-volume gap is assumed as
input, and no claim of Lean verification, open-boundary coverage, or arbitrary aspect ratios
is made.

## Layout

- `src/SquareIsing.jl`: exact row models, transverse symmetries, and RP factors.
- `src/CubeRecurrence.jl`: exact three-component recurrence and interpolation-window bounds.
- `src/Certificate.jl`, `certify.jl`, `check.jl`: certificate generation and replay.
- `scripts/build_spectra.jl`, `scripts/scan.jl`: verified spectra and numerical seed search.
- `test/runtests.jl`: independent geometry, arithmetic, spectral-moment, and rejection tests.
- `METHOD.md`, `IDEAS.md`: proof and route analysis.
- `results/`: relative symlink to `../../../data/models/square_tfim/`, containing certificates,
  spectrum manifests, compressed verified spectra, and numerical scans.

Model-local `src/`, `scripts/`, `certify.jl`, and `check.jl` hashes are enforced. Core hashes
are provenance only; the checker reruns core arithmetic. Spectrum files are hash-pinned by
the certificates. All can be regenerated in a few minutes; no immutable source exceptions
are needed.
