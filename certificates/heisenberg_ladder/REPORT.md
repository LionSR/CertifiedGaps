# Spin-1/2 two-leg Heisenberg ladder

**Certified:** at J_perp=4, Δ_L ≥ 1.102403170532509 for every even L ≥ 4 rungs
(2L ≥ 8 spin sites, with the site count divisible by four). The exact binary
lower bound is the `.bound` value in `results/certificate.json`.

## Model, scope and conventions

H_L = sum_{j=1}^L sum_{a=1}^2 S_{a,j}.S_{a,j+1}
+ J_perp sum_{j=1}^L S_{1,j}.S_{2,j}, with periodic j and J_parallel=1.
One cell is one rung, of dimension 4. At L=2 the periodic leg bonds are counted
twice, exactly as in the core's `ring_bonds` convention. The scan uses J_perp =
4, 7/2, 3, 5/2, 2, 3/2; all are exact rational inputs to matrix construction.

**The seed must have even n0.** The source scope note (Appendices D–F) explicitly
requires even powers of the self-adjoint spatial transfer operator. Therefore
rings up to 10 rungs allow valid n0=2,4 only. We compute n0=3,5 as diagnostics,
but neither can certify this one-rung-cell model. Odd rings are frustrated and
the checkerboard rotation is not periodic there. In particular a negative odd-n
spatial defect or a spurious odd seed cannot define the requested threshold.

The threshold below is a **numerical finite-search threshold**, not a physical
critical coupling, not proof of monotonicity in J_perp, and not proof that a
failed parameter has no seed at larger sizes or off the beta grid.

## Exact reflection positivity

In the spin-z basis put

```
X = (1/2)[0 1; 1 0], A = i Sy = (1/2)[0 1; -1 0],
Z = (1/2)[1 0; 0 -1], u = 2A = exp(i pi Sy).
```

All entries are rational. The original two-spin bond is
h = X⊗X - A⊗A + Z⊗Z. Rotation of exactly one endpoint by u gives
h' = -X⊗X - A⊗A - Z⊗Z = -Sx⊗Sx + Sy⊗Sy - Sz⊗Sz.
The stated checkerboard rotation rotates one endpoint of every leg and rung
bond for even L. On alternating rungs the rotations are u⊗I and I⊗u.
Both give the same on-rung term D=J_perp h', which is real symmetric.
Both orientations also give the same inter-rung bond

```
-h_bond = sum_{C in {X⊗I,A⊗I,Z⊗I,I⊗X,I⊗A,I⊗Z}} C⊗C.
```

Every coefficient is exactly 1, hence nonnegative. `rp.jl` checks both
orientations, the original-to-rotated identity, the sum, and symmetry of D
using Rational{BigInt}, for every requested coupling. Real antisymmetric A is
allowed; Hermiticity of C is not required by the stated RP hypothesis.
The rotated one-rung Hamiltonian is translation invariant. This establishes
the real self-adjoint transfer-operator hypothesis of the supplied criterion.

Ground-state uniqueness also has a direct finite-volume argument. For positive
couplings the connected bipartite even ladder has a Marshall sign transform
with negative exchange off-diagonals. In each fixed magnetization the spin-swap
graph is connected, so Perron–Frobenius makes its lowest vector strictly positive
and simple. In total magnetization zero, a product of rung singlets has, after
the same sign transform and an overall phase, nonnegative nonzero amplitudes.
Its overlap with the lowest vector is positive. Since SU(2) symmetry and
simplicity make that lowest vector a definite-total-spin vector, the overlap
with a singlet forces total spin zero. Any other global ground state would
have a magnetization-zero member of its SU(2) multiplet, contradicting this
simplicity. This proves uniqueness for each even ring at the stated positive
couplings; a uniform gap is established only by the verified seed below.

## Full dense spectra and exact counting

The ED is in the original unrotated spin basis. An exact unitary rotation
preserves the partition function. Bits encode the 2L spins. Exact rational
matrix elements are diagonal ±J/4 and spin exchange J/2.

The abelian symmetry group is rung translation C_L, leg exchange C_2, and,
only at zero magnetization, global spin reversal C_2. Orbit representatives
and stabilizers are computed with integer bit permutations. A character with
labels (k,p,q) is retained precisely when it is trivial on every stabilizer;
this test is modular integer arithmetic, with exponent 2kt+L(pa+qb) mod 2L.
Normalized orbit sums are an orthonormal basis of each character subspace.
Their matrix elements are sums of exact rational transitions times
sqrt(orbit_size_source/orbit_size_target) times a root of unity.

Opposite magnetizations are isospectral under global spin reversal; k and L-k
are isospectral under complex conjugation of the real original Hamiltonian.
Their exact multiplicities are 2, except at fixed sectors. Every retained
block is diagonalized **densely**. No partial spectrum or Lanczos estimate is
used. For each L, the code asserts sum(block_dimension × multiplicity)=4^L,
and also asserts the orbit sizes sum to binomial(2L,m) in each charge sector.
Each JSON spectrum records every block label, dimension, multiplicity and time.
Small-L smoke tests compare all eigenvalues against an independent full matrix
and the core's separate `Chain`/`sector_hamiltonian` implementation, through L=5.
They also check the analytical normalized second trace moment.

For numerical scans, rung reflection times complex conjugation is an
antiunitary involution preserving each momentum sector and commuting with H.
Its fixed real basis converts the complex block into a real symmetric block
of the same dimension. `realify!` constructs that unitary orbit-pair change of
basis and checks the imaginary remainder. The last two scans use
`sparse_scan_basis.jl` to carry out this transformation before densifying;
`test_lowmem.jl` repeats the same independent full-spectrum tests. These
optimizations affect only numerical scans; the certified path encloses the
original symmetry blocks directly, including complex momenta.

## Implemented enclosure theorem

The rigorous run evaluates every matrix entry using IntervalArithmetic:
rational coefficients, square roots of rational orbit-size ratios, and sinpi/
cospi of rational angles are all enclosed outward. It selects a Hermitian
floating matrix Hhat from the upper triangle; exact Hermiticity fixes the
lower triangle. The maximum row sum of entrywise complex errors, using
|Re error|+|Im error|, gives h_error >= ||H-Hhat||_2.

Let Q, Lambda be the complete floating dense eigendecomposition. The actual
residual R=Hhat Q-Q Lambda and orthogonality defect E=Q*Q-I are bounded as
follows. The BLAS products are evaluated in Float64, then rigorous interval
Frobenius norms and explicit roundoff allowances enclose their errors. With
u=2^-53 and gamma=(8d+32)u/[1-(8d+32)u], the implemented upper bounds are

```
r   = ||R_fl||_F + gamma (||Hhat||_F ||Q||_F + ||Q||_F ||Lambda||_F) + 2^-1000
eta = ||E_fl||_F + gamma (||Q||_F^2 + sqrt(d)) + 2^-1000.
```

The factor 8d+32 conservatively covers real operations in a length-d complex
dot product, diagonal multiplication and subtraction. The standard
IEEE round-to-nearest dot-product forward-error bound is componentwise
|fl(AB)-AB| <= gamma |A||B|; Cauchy–Schwarz gives the displayed Frobenius
bounds. Every norm and scalar bound is accumulated with outward interval
arithmetic, not accepted from an unverified floating norm. The additive
2^-1000 dominates subnormal absolute errors at these block dimensions.
There is no overflow. OpenBLAS uses ordinary dot-product arithmetic, with
BLAS thread count fixed at one. Thus higher precision cubic products are
unnecessary: this is the rigorous-norm-upper-bound alternative in rule 3(c).

Here is the Weyl/polar-factor theorem and its proof, to make the enclosure
independent of an unstated eigensolver accuracy assumption. If eta<1, the
polar decomposition Q=US has unitary U and

```
||Q-U||_2 <= eta/(1+sqrt(1-eta)) = d_Q.
||Hhat U-U Lambda||_2 <= r + (||Hhat||_F+max|lambda_i|) d_Q.
```

The first inequality follows by applying |sqrt(x)-1| to the eigenvalues of
Q*Q in [1-eta,1+eta]. The second follows by adding and subtracting Hhat Q
and Q Lambda. Unitary conjugation and the Hermitian Weyl inequality then
give, for every sorted eigenvalue, **with multiplicity**,

```
|lambda_i(H)-lambda_i| <= rho
rho = h_error + r + (||Hhat||_F+max|lambda_i|) d_Q.
```

`verified_eigen` implements precisely these outward operations. Each block
stores rational upper bounds for residual, orthogonality defect, matrix error
and rho. All expanded eigenvalues and radii are stored as exact rational
strings in `verified_L*.json` (the exact dyadic numbers denoted by Float64).
This establishes a complete enclosure without trusting LAPACK's error estimates.

For each beta, interval exponentials of every eigenvalue interval are summed,
with an exact floating shift for stability, then logarithms give [lo,hi] for
log Z. Consequently

```
P_s >= exp(lo(logZ_2n_beta) - 2 hi(logZ_n_beta))
P_t >= exp(lo(logZ_2n_2beta) - 2 hi(logZ_2n_beta)).
```

The implemented interval `-expm1` expressions return upper defect endpoints;
these are rounded up again to rationals on the 10^-12 grid. Those rationals
and exact rational beta are passed to `certified_gap_bound`, with its default
J=48 outward-rounded recurrence and tail correction.

The checker repeats exact RP, sector dimension accounting, source hashes,
and the bound from stored rational seeds. It includes `recheck_seeds.jl`, which
recomputes log partition functions and seeds from all stored eigenvalue
enclosures. It then performs fresh complete verified 4/8-rung eigensolves,
recomputes every residual and orthogonality allowance, and checks that the
stored rational seeds still dominate the newly computed defects. Thus the
checker does not merely trust saved norm bounds. `validate.jl` asserts that
this full checker completes in under two minutes. `certify.jl` regenerates
all certificate JSON files.

## Environment and timing scope

Julia 1.12.4 and the checked-in dependency manifest were instantiated and
precompiled successfully. The approved juliaup endpoint returned HTTP 403;
the official Julia 1.12.4 binary download succeeded. Julia and its depot live
under `/tmp`, while every task script and result lives in this certificate
directory. All runs set OPENBLAS_NUM_THREADS=1 and BLAS.set_num_threads(1).
Parameter processes share the cloud environment's four CPU cores. No work was
moved to the user's laptop. A reported connection interruption was checked;
the existing Julia processes and saved certificate remained available.

The scan timings below are for completed driver invocations, which reuse
saved smaller-ring spectra. Per-ring JSON files preserve their original
compute times. Setup, compilation warmup, and interrupted exploratory complex
10-rung runs are not included in those scan times. Those exploratory runs were
stopped before restarting with a real reflection basis; no duplicate active
scan of the same coupling was left running. One later dense-intermediate
J=7/2 10-rung attempt was killed by the memory limit (exit 137, confirmed
cgroup OOM event). Its small-ring checkpoints and the completed certificate
survived. Only the missing large-ring calculation was queued for retry with
the tested sparse-to-real construction. Final results below use completed
runs, never partial spectra from the killed process.

## Results

Results and executed validation timings are appended below by `summarize.jl`.

## Current-main reproduction

The original cloud scans used repository base `ef0a3e6`. Upstream `22b2759`
renamed the unchanged mathematical core to `CertifiedGaps`. Publication uses
that current base in an isolated worktree, changes only this certificate
directory, and reruns the complete rigorous certificate and checker there.
Historical numerical scan metadata retains its original source hashes and
package name. These scans are numerical evidence; they are not promoted to
certificates by the namespace update. The requested `certificates/` location
is retained even though the newer repository guide prefers `models/`.

## Separate current-head audit

`audit.jl` passed in 12.763 seconds. It did not read stored radii:
it replayed all 256 and 65,536 eigenvalues from fresh verified dense solves,
checked the exact normalized trace and second trace moment of the ladder,
and independently recomputed the three log partition functions with 256-bit
MPFR directed rounding. Both resulting defect upper bounds fit within the
stored rational seeds, and the recurrence reproduced the same bound.

This is independent partition arithmetic and exact trace-moment checking,
not an alternative implementation of every symmetry block or norm theorem.
Fresh eigensolves share the symmetry builder and residual-bound implementation.
The abstract repository theorem and ordinary IEEE/OpenBLAS arithmetic model
remain assumptions; no rigorous thermodynamic reference gap or threshold is
claimed. The pre-rename validated snapshot is preserved byte-for-byte under
`legacy/`, including its matching certificate source files.

The migration audit compared all 25 core source files at `ef0a3e6`
and `22b2759`: after replacing the package name, their bytes match exactly.
The recurrence implementation itself is byte-identical.

<!-- GENERATED RESULTS -->

### Numerical scan status (draft)

The even-seed grid is complete through n0=4. Full 10-rung scans are complete
for J=4,3,5/2,2,3/2; only J=7/2 is still running. No unfinished 10-rung
result is included in a completed-scan claim.

| J | Best (n0,beta) | Numerical recurrence bound |
|---:|:---|---:|
| 4.0 | (4, 2) | 1.1024044365 |
| 3.5 | (4, 2.25) | 0.7786489922 |
| 3.0 | (4, 2.5) | 0.4541898038 |
| 2.5 | (4, 3) | 0.0424345905 |
| 2.0 | none | none |
| 1.5 | none | none |

The numerical finite-grid threshold lies between 2.4881973267 and 2.4882354736.
This is not a certified coupling threshold or a physical critical point.
