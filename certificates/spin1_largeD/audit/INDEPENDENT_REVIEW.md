# Independent mathematical and implementation review

Reviewed certificate implementation at commit
`4356117476ec92fa0ff9fd002f5259f9851cf94b`, with current core revision
`22b2759e0790197b83e2571b465aa8dc69ebc0be` (`CertifiedGaps`). The reviewer was
delegated separately from the worker performing the fresh dense spectral replay.
No original results, archived sources, or shared core were changed by this review.

The additional reviewer script passed **129 exact Cartesian/Sz checks** and
**90 residual/orthogonality diagnostic checks** in 4.5 and 5.4 seconds of test
execution respectively (excluding Julia startup/compilation). See
`reviewer-checks.log`. Reproduce from repository root with Julia 1.12.4:

```
OPENBLAS_NUM_THREADS=1 julia --project=julia certificates/spin1_largeD/audit/reviewer_checks.jl
```

On this environment Julia is `/tmp/julia-1.12.4/bin/julia` and the installed
dependencies require `JULIA_DEPOT_PATH=/tmp/spin1-julia-depot`.

**Verdict:** no mathematical or implementation defect was found in the reviewed
model-specific enclosure pipeline. The two uniform conclusions are supported
under the repository's two-purity theorem and the stated IEEE/BLAS arithmetic
model, subject to successful fresh dense replay recorded separately. This review
does not itself assert completion of that expensive replay. D=1603/500 is the
smallest successful sampled value, not a global threshold.

## Exact model and RP normalization

The supplied integer Cartesian matrices give Hermitian spins `Sa=-i Ga`,
`[Sx,Sy]=i Sz`, and `sum Sa^2=2I`. Their two-site sum equals
`P-|Omega><Omega|`, with the **unnormalized** vector
`Omega=|xx>+|yy>+|zz>`. Consequently `-h=sum Ga⊗Ga` has coefficients exactly
one. Each Ga is real, so its second copy is also its complex conjugate. The
onsite matrix is exactly `D diag(1,1,0)`. There is no missing factor of two or
three in this RP identity.

The numerical Hamiltonian is built in the Sz basis. The independent exact check
in `reviewer_checks.jl` now bridges the two bases rather than comparing two
uses of `transitions()`. Its unnormalized Sz eigenvectors are the columns of

```
V = [ i  0 -i
      1  0  1
      0  i  0 ],       V* V = diag(2,1,2).
```

Here `V*` means adjoint. Dividing the columns by their norms gives m=-1,0,+1;
the middle-column phase makes the ladder amplitudes positive. Transforming the
Cartesian bond with integer-complex arithmetic and exact normalization gives
diagonal `m1*m2` and every permitted opposite ladder move with amplitude one.
The non-square normalization cases have identically zero numerators. The test
also checks the onsite transformation and the L=2 convention of two periodic
bonds, for both certified D values.

## Sectors and block construction

Charge conservation, translation, and spin inversion are exact symmetries of
the rational product-basis transitions. Charge inversion pairs q with -q.
In q=0, the commuting group is translation times spin inversion. An orbit is
present in a character sector precisely when its character is trivial on both
the translation stabilizer and, where applicable, the stabilizer combining
translation with spin inversion. `active_orbits` implements these conditions
with integer modular arithmetic.

Distinct orbit sums have disjoint support. A normalized projected orbit vector
has coefficient magnitude `1/sqrt(size(a))` on its support. Equivariance of H
then gives the representative transition formula

```
h * sqrt(size(a)/size(b)) * exp(-2π i k u/L) * parity^f.
```

This is the factor used in `block`; interchanging the size ratio would be wrong.
The conjugate momentum sectors have equal spectra, because H is real and spin
inversion parity is real. The charge and momentum multiplicities therefore
give a direct-sum decomposition, not merely a dimension coincidence. The
weighted dimension assertion `sum(mult*dim)==3^L` is an additional check.
Every block in the rigorous branch is passed to a dense Hermitian eigensolver;
the partial Krylov branch is separate and never used there. The L=12 recorded
weighted dimension is 531441 and the largest dense block is 5798.

Exact rational transitions, phases (including an interval enclosure of π), and
square roots are independently accumulated as intervals. Using the independent
upper triangle and its adjoint is justified by the exact Hermiticity of the
projected operator. The entry bound `|Re(error)|+|Im(error)|` dominates its complex
modulus. Summing it into both affected rows gives a maximum row-sum upper bound
`c` on the spectral norm of the Hermitian construction error. The physical
norm bound `L(2+|D|)` is valid because the spin-one bond has eigenvalues -2,-1,1
and the onsite projector has norm one. Thus `M=L(2+|D|)+c` bounds the constructed
matrix A.

## Residual, orthogonality, and eigenvalue radii

The eigensolver outputs are not assumed exact. The code computes full residual
and orthogonality matrices and enlarges their Frobenius norms. Under the usual
binary64 arithmetic model, a complex dot product of m terms is covered by
`gamma_(8m+32)`; this also dominates the extra real diagonal multiplication and
final subtraction. Componentwise matrix-product bounds imply

```
||fl(AQ)-AQ||F <= gamma ||A||F ||Q||F
||fl(Q*Q)-Q*Q||F <= gamma ||Q||F^2.
```

Since `||A||F <= sqrt(m) M`, these are precisely the terms used, together with
the `Q Lambda` and identity/subtraction terms. The absolute-value matrix product
bound is legitimate here: its Frobenius norm is at most the product of the two
input Frobenius norms. No unsupported spectral-norm bound on an entrywise
absolute-value matrix is being used.

For N entries, `sum(abs2,X)` is nonnegative. Its real operations are bounded
conservatively by `4N+16`; division by `1-gamma_(4N+16)` converts the computed
sum into an upper bound. Taking an outward-rounded square root gives a valid
Frobenius upper bound. At m<=5798 both gamma denominators are positive. The
additions `1e-290` before the square root and `1e-280` in the residual bounds
dominate accumulated subnormal errors at these dimensions: even a crude
10^13 elementary-operation count times the binary64 normal/subnormal boundary
is below 3e-295 (the squared-norm sum itself uses fewer than 2e8 operations,
whose corresponding bound is below 5e-300). No unverified higher-precision
eigensolver accuracy is needed.
The supplementary 256-bit residual tests are numerical diagnostics of these
inequalities, not the source of their rigor.

For `e<1`, the polar factor U of Q is unitary and
`||Q-U|| <= e/(1+sqrt(1-e))`. Consequently the code's

```
eta = c + r + (M+max(abs(lambda))) * e/(1+sqrt(1-e))
```

bounds `||Hblock-U Lambda U*||`. Both matrices are Hermitian, so min-max/Weyl
encloses every sorted eigenvalue, with multiplicity. This also handles exact
degeneracies; it is not an eigenvector-by-eigenvector matching argument.

## Partition direction, recurrence scope, and uniqueness

An eigenvalue interval is `[lambda-eta,lambda+eta]`. For positive beta its lower
endpoint gives an upper exponential weight and vice versa. Summing all weights,
taking an interval log, and undoing the same exact binary energy shift encloses
each log Z. The **lower** log-ratio endpoint gives the lower purity bound;
therefore `sup(1-exp(logZ2-2logZ1))` is the required **upper** defect bound.
The conversion of that Float64 endpoint to `Rational{BigInt}` is exact.
Both spatial and thermal ratios in `seed.jl` have the correct direction and
normalization. Upper log-Z bounds alone would not suffice; these are two-sided
intervals.

For clarity, the all-even-length scope does not follow just by doubling the
seed lengths. RP supplies real spatial-transfer eigenvalues; even moments are
nonnegative. Write x=L/n in [1,2]. Moment log-convexity at beta and monotonicity
of the corresponding ell-p norms at 2beta give

```
log Ptau(L,beta) >= (x/2) log Ptau(2n,beta)
                       + (2-x) log Ps(n,beta)
                 >= log Ptau(2n,beta) + 2 log Ps(n,beta).
```

Thus every even L in [n,2n] has thermal defect at most
`z(s,t)=1-(1-t)(1-s)^2`. The consecutive dyadic bands starting at even n=6
cover all even L>=6 sites. Odd n=7 (L=14 spectrum) would not directly be an
eligible seed. Coupling and squaring give the asymmetric recurrence implemented
in the core. Entering the admissible square after some iterations is valid;
initial defects need not individually be below 0.12731.

The core checks z<1/2 at every finite step and subtracts its rigorous tail
allowance. Monotonicity before entrance to the square can also be seen directly:
with r=z/(1-z)<1, the next spatial defect is r^2/2 and hence the next combined
defect is at least r^2-r^4/4 >= r^2/(1+r^2). This makes the gap rates
nonincreasing. The final tail bound is the unchanged shared-core theorem; this
review is not a fresh formal verification of the core or its dependencies.

Thermal purity greater than 1/2 implies the largest Gibbs weight exceeds 1/2,
so the ground state is unique. Moreover the ratio of the next Gibbs weight to
the largest is at most z/(1-z), giving the gap estimate. Uniform uniqueness
therefore follows from the moment/recurrence argument; it is not inferred by
extrapolating the two finite-ring gaps. The strictly positive finite-ring gap
intervals in `make_seed` are additional checks only.

## Replay distinction and remaining audit limits

The historical and fast/current checkers verify fingerprints and recompute
partition intervals and recurrence bounds **from stored radii**. They do not
rebuild dense matrices, recompute eigenvectors, or recompute residual and
orthogonality bounds. Therefore those checkers alone are not the independent
spectral replay requested in the final audit. The parent worker's separate
replay must rebuild both L=6 and L=12 spectra at D=4 and D=1603/500 and retain
new enclosure evidence without overwriting the originals.

The code/proof audit found no required mathematical correction. Acceptance
still depends on the separate replay outcome. The trusted base remains Julia,
IntervalArithmetic's outward rounding, the standard IEEE/BLAS operation model,
and the repository's moment/recurrence theorem. The BLAS executable and the
entire theorem have not been formally verified by this review. Search failures,
the ED reference gap, ratios, and any claim concerning all smaller real D are
not certified by the two successful seeds.

## Fresh replay driver and assembly addendum

The separately reviewed `replay_worker.jl` (SHA-256
`3ed2250c52ba9ba4f1d01f987f0a2a70bd84a53890447262376351d5e2c76d45`)
reads no historical spectrum or radius file. It starts from exact rational D,
reconstructs all matrix blocks with interval construction errors, and computes
new full eigenvectors, residuals, orthogonality defects, and radii. Its two
charge groups `{1,4,...,12}` and `{0,2,3}` are disjoint and exhaust q=0,...,12.
The mixed group additionally rebuilds the entire L=6 spectrum. No dense jobs
were started by the independent reviewer.

The assembler enumerates sectors independently from orbit data and requires
exactly 92 nonempty L=12 blocks, weighted dimension 531441, and 29 L=6 blocks,
weighted dimension 729. It verifies sector identities, multiplicities, full
spectrum lengths, finite numbers, orthogonality defects below one, D values,
source hashes, core tree, and disjoint complete worker manifests. Added
hardening checks also reject missing radii, incorrect L=6 sector metadata,
incomplete worker listings, dirty tracked core files, and a core tree differing
from the audited current revision. The resulting record hashes the actual
current core files, not just its Git tree name.

Fresh log-partition intervals, rational upper defects, finite-ring gap
intervals, and a positive recurrence bound are recomputed **before** reading
historical results. Only then are the historical spectra and certificate read
for comparisons. The published nine-decimal lower bounds are supplied as exact
rationals and compared against the exact binary value of the new bound. The
historical center/radius comparisons are diagnostics and do not supply a
certificate premise. `assemble(...; save=false)` performs the same validation
and returns the fresh record without writing a result, enabling a later checker.

The assembly implementation was statically reviewed and syntax-checked while
the worker jobs were still running. Completion and actual replay outputs must
be read from the separate final replay record; this addendum does not claim
those computations had already finished at review time.

## Final independent verdict after completed replay

The reviewer subsequently inspected both final `audit/results/certificate_D*.json`
records, `assembly.log`, `fresh-check.log`, the worker completion logs and
manifests, and the final assembler/checker sources. An independent hash pass
rechecked **252 file-hash entries**: for each certificate, 95 fresh input files,
four source files, and 27 current-core files. Every hash matched. The frozen
worker hash above remained unchanged; the final assembler hash was
`845c79469447739408a4b59bec7f7dd32af4399a0f2d6bb6ed79f5e47163baa7`.

All four fresh worker jobs completed in 1265.31–1489.78 seconds (under 25
minutes each). At each D the heavy worker supplied 64 blocks with weighted
dimension 254592, and the mixed worker supplied 28 L=12 blocks with weighted
dimension 276849, giving the required 92 blocks and dimension 531441. The
mixed workers also rebuilt L=6. Both records match the audited current core
tree `5ebfb34b3f28f37a6fe600857b387dfe887de8b8`.

The fresh computations reproduced all historical L=6 and L=12 eigenvalue
centers and enclosure radii bit for bit, both exact rational seed pairs, and
both raw recurrence outputs:

| D | beta | Reproduced lower-bound endpoint | Safe displayed bound |
|---|---|---:|---:|
| 4 | 9/4 | 0.6946941929576096 | 0.694694192 |
| 1603/500 | 381/125 | 0.0022134568960786357 | 0.002213456 |

The fresh checker passed **240/240 checks**, with logged wall time 17.896
seconds. Its role remains validation of the newly stored evidence; the dense
worker jobs, not this fast checker, recomputed spectral radii from model
inputs. The documented initial driver launch failures are separate from these
successful completed jobs and do not supply certificate inputs.

**Final verdict: accepted within the stated trusted base.** No outstanding
model-specific mathematical, enclosure, completeness, provenance, or replay
issue was found. Both bounds apply to every even L>=6 one-site cells (sites),
and uniqueness follows as explained above. This resolves the replay condition
in the earlier provisional verdict without changing that historical record.
The trusted-base limits remain the repository's moment/recurrence theorem,
Julia and interval rounding, and the standard IEEE/BLAS error model; the
review does not formally verify those foundations or provide a second
independently implemented large-matrix eigensolver. D=3.206 is still only the
smallest successful sampled D. Search failures, physical gap estimates and
ratios, and a global threshold claim remain numerical or unproved.
