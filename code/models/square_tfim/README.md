# Square-lattice transverse-field Ising model

The spacetime bootstrap certifies a unique ground state and a uniform gap on every even
square torus of side at least four. Both spatial periods grow. A direct-temperature seed
prefix improves the use of complete verified spectra of at most sixteen spins; the
subsequent three-direction recurrence and spatial interpolation are unchanged.

| parameter $g/J$ | seed $(n_0=m_0,\beta_0J)$ | certified gap $\ge$ | lengths covered |
|---:|---|---:|---|
| 15 | $(4,47/100)$ | $0.74208J$ | even square side $L\ge4$ |
| 16 | $(4,141/320)$ | $1.72356J$ | even square side $L\ge4$ |
| 18 | $(4,143/360)$ | $3.39128J$ | even square side $L\ge4$ |
| 20 | $(4,29/80)$ | $4.97710J$ | even square side $L\ge4$ |
| 24 | $(4,5/16)$ | $8.14192J$ | even square side $L\ge4$ |

## Model

With Pauli matrices and periodic indices,

```math
H_{n,m}=-J\sum_{x,y}(Z_{x,y}Z_{x+1,y}+Z_{x,y}Z_{x,y+1})-g\sum_{x,y}X_{x,y}.
```

The computation sets $J=1$. Lengths count individual lattice spacings; an $L$ by $L$
torus has $L^2$ spins. At period two each bond in that direction is counted twice, as
required by the periodic sum and transfer trace. For $L\ge4$, each undirected bond in
the square Hamiltonian occurs once. The new bounds do not claim coverage of $L=2$.

## Reflection positivity

In the physical basis a crossing bond has $-h=Z\otimes Z$. A Hadamard rotation on every
spin gives the complete-row interaction $-\sum_y X_y\otimes X_y$ and real symmetric
on-site term $-\sum_y X_yX_{y+1}-g\sum_y Z_y$. `SquareIsing.row_rp` and the core
`verify_rp` check the integer factors $C_y=X_y$ exactly.

Finite Trotter transfers in both spatial directions have the form $D^{1/2}KD^{1/2}$,
with real symmetric positive K. Their squaring and moment interpolation inequalities
pass to the exact quantum partition function. The proof is in [METHOD.md](METHOD.md);
no numerical Trotter or tensor truncation is used.

## Seeds and certified bounds

At half the effective seed inverse temperature, $\beta=\beta_0/2$, the verified spectra directly
enclose the defects of

```math
u=P_y(4,2,2\beta)=\frac{Z(4,4,2\beta)}{Z(4,2,2\beta)^2},\qquad
 t=P_\tau(4,4,\beta)=\frac{Z(4,4,2\beta)}{Z(4,4,\beta)^2}.
```

The prefix uses these two bounds to obtain the staggered defects of
$P_x(4,4,\beta_0)$, $P_y(8,4,\beta_0)$, and $P_\tau(8,8,\beta_0)$:

| $g/J$ | $\varepsilon_x$ | $\varepsilon_y$ | $\varepsilon_\tau$ |
|---:|---:|---:|---:|
| 15 | 0.008988535251 | 0.044342296190 | 0.020778805396 |
| 16 | 0.006829707080 | 0.032764009742 | 0.016668809184 |
| 18 | 0.004282379581 | 0.019758812120 | 0.010617692121 |
| 20 | 0.002839873795 | 0.012752261131 | 0.007384511624 |
| 24 | 0.001433660730 | 0.006222333403 | 0.003739202719 |

These are exact terminating-decimal upper bounds. The face/cube identities and the usual
squaring lemma prove the prefix; no new transfer assumption or larger spectrum is needed.
The complete 4 by 4 spectrum has 65,536 levels with exact multiplicities, and its largest
symmetry block has dimension 2392. The 2 by 2 spectra remain in the common spectral record
for comparison with the earlier certificates; the new prefix needs only 4 by 2 and 4 by 4.

Twelve exact three-component updates precede the quadratic tail. The stronger bounds use
a rational grid of denominator $10^{8000}$ so that rounding does not set the tail's scale.
Every earlier interpolation window is included, giving all even sides rather than only
doubled sizes. The one-dimensional `certified_gap_bound` is not used.

Each certificate also proves a closed field interval using a finite-volume Weyl estimate
for $H/g$. This is a uniform parameter bound, not interpolation of numerical eigenvalues:

| interval of $g/J$ | uniform gap $\ge$ |
|---|---:|
| $[14.999,15.001]$ | $0.51384J$ |
| $[15.999,16.001]$ | $1.55379J$ |
| $[17.999,18.001]$ | $3.23759J$ |
| $[19.999,20.001]$ | $4.82167J$ |
| $[23.999,24.001]$ | $7.97524J$ |

The `legacy` section preserves the original and refined certificates at $g/J=18,20,24$,
including their numerical values. Their former best bounds were $0.74378J$, $2.51376J$,
and $5.52173J$. The current primary scope is $L\ge4$, even though the legacy certificates
also cover the repeated-bond $L=2$ torus. No continuous bootstrap regime $g/J\ge15$ is claimed.

The separate ferromagnetic Ising graph argument gives
$`\Delta\ge2\sqrt{g(g-4J)}`$ for $g>4J$, and is substantially stronger numerically.
The results here improve the finite-data spacetime bootstrap, not the known high-field
physical regime. [METHOD.md](METHOD.md) gives the comparison and the archived proof reference.

## Reproduce

From the repository root, with Julia 1.12 and the core environment:

```sh
export OPENBLAS_NUM_THREADS=1
julia --project=code/julia code/models/square_tfim/certify.jl
julia --project=code/julia code/models/square_tfim/check.jl
julia --project=code/julia code/models/square_tfim/check.jl --full
julia --project=code/julia code/models/square_tfim/test/runtests.jl
```

The default drivers process all five fields. Supply one field, for example `check.jl 15 --full`.
`certify.jl --reuse` reuses model-checked, hash-pinned spectra. Without it, every verified
spectral block is regenerated. The full replay regenerates spectra independently of the
stored enclosures, then verifies that their fresh bounds imply the saved seed limits.
On the development Apple Silicon machine with one BLAS thread, the five-point fast
replay takes 53.8 seconds, reused-spectrum regeneration 53.6 seconds, and the 1,161 model
tests 176.4 seconds in their test sets. These timings exclude process startup and
precompilation. Full spectral replay timings and detailed checks are in [METHOD.md](METHOD.md).

For floating-point temperature searches only:

```sh
julia --project=code/julia code/models/square_tfim/scripts/refine_scan.jl --direct
julia --project=code/julia code/models/square_tfim/scripts/refine_scan.jl
julia --project=code/julia code/models/square_tfim/scripts/moderate_fields.jl
```

The first command scans the direct prefix; the second reproduces the earlier three-point
refinement. Neither scan is a certificate. The general Lean theorem can be checked with
`cd code/lean && lake build CertifiedGaps.CubeChain`.

## What is rigorous

The spatial transfer proof, exact RP identity, full spectral enclosures, interval partition
sums, direct-temperature prefix, parameter perturbation estimate, and exact rational cube
chain establish the bounds together. The spectral trusted base is IEEE binary64
round-to-nearest with gradual underflow, the a priori dot-product rounding bound for BLAS
(including arbitrary summation order and FMA), verified residual/orthogonality estimates,
and `IntervalArithmetic`. Final logarithms use the exact rational series in `ExactArithmetic`.

The default checker verifies model-local source hashes, exact row models and RP factors,
spectral hashes and multiplicities, fresh partition enclosures, both new chains, every
parameter-interval quantity, finite-volume energy enclosures, and the legacy records.
The tests include independent geometry and spectral moments, an independent audit of the
chain inequalities, and deliberate changes to seeds, gaps, field intervals, scope, and hashes.
Long certificate rationals are already quoted exact strings. The model reader uses the
core's typed JSON reader rather than its regex-based integer-token conversion, and rejects
unquoted scientific numbers.

Lean's `CertifiedGaps.Spacetime.gapGT_of_cube_chain` proves the general criterion from
explicit scalar moment, spatial interpolation, and seed hypotheses. The Ising transfer
construction, new finite prefix, parameter perturbation proof, and numerical seeds remain
outside Lean; these Julia chains have not been imported as Lean witnesses. This is not an
end-to-end Lean proof of the Ising Hamiltonian's gap. The field scans remain floating-point
diagnostics. Ordered-phase splitting and the gap above the lowest pair are different
quantities; the current unique-ground-state criterion does not cover that phase.

## Layout

- `src/SquareIsing.jl`: exact row models, transverse symmetries, and RP factors.
- `src/CubeRecurrence.jl`: direct-temperature prefix, exact recurrence, and window bounds.
- `src/Certificate.jl`, `certify.jl`, `check.jl`: one generation/replay path for all points,
  closed parameter intervals, and preserved legacy records.
- `scripts/`: spectral generation and clearly labeled numerical diagnostics.
- `test/runtests.jl`: independent geometry, arithmetic, moment, and rejection tests.
- `METHOD.md`, `IDEAS.md`: proofs, comparisons, route analysis, and earlier diagnostics.
- `results/`: relative symlink to `../../../data/models/square_tfim/`.

Model-local `src/`, `scripts/`, `certify.jl`, and `check.jl` hashes are enforced. Core hashes
are provenance; the checker reruns the arithmetic. No immutable source exceptions are needed.
