# Independent acceptance review

Review date: 2026-10-11. Scope: the three stored general-reflected-RP
certificates, their mathematical justification, and regeneration of their
spectral witnesses. This review does not extend the parameter or size grid.

## Verdict

**PASS: all three general internally reflected RP certificates.** The final
independent acceptance run passed **4012/4012 checks in 145.5 seconds**,
regenerating all six spectra and rigorously validating containment of every
original stored eigenvalue interval. **FAIL: the original stronger literal
real `C tensor C` hypothesis.** That condition is not satisfied and must not
be claimed. There is no remaining acceptance blocker for the three general-RP
certificates under the explicitly stated classical binary64 OpenBLAS trusted
arithmetic model. The result concerns the specified physical Hamiltonian
without rescaling energies or changing boundary conditions.

| delta | accepted stored uniform bound | independently regenerated refined bound |
|---|---|---|
| 1/2 | 0.5326577398571627 | 0.5347797088100107 |
| 3/10 | 0.17291646480187198 | 0.17320162195431488 |
| 1/4 | 0.057966591686137485 | 0.058229681460245694 |

The refined values are acceptance evidence; the production certificates and
their reported bounds are unchanged. Both columns use the same general-RP
theorem and the range every even `L>=4` cells.

## Exact theorem chain and normalization

All repository paths below are relative to the repository root; line references
refer to the source snapshot whose hashes are checked by the stored checker.

- `outputs/thermal-mpo-gap-scope.tex:75` supplies locality and translation
  invariance, applied here to identical two-spin cells of dimension four.
- The same file, lines 85–106, explicitly permits an antiunitary
  `Theta=V K`, `Theta^2=+/-1`, the bond form
  `-h=sum c C tensor Theta(C)`, and `Theta(D)=D`. It does not require real
  identical factors. `rp.jl:26–39` verifies the actual rotated physical bond,
  real channels, involution, and onsite invariance in exact rational arithmetic.
- The mirror lemma is at scope lines 130–154. Its proof at lines 148–153
  preserves chronological order: antiunitary conjugation is an antilinear
  algebra automorphism, not a transpose or an adjoint. For `V=F_cell`, it
  exchanges the left/right endpoint insertions while fixing `Drot`. All
  product-basis insertions are real, so this particular kernel is real
  symmetric as well as Hermitian. There are no fermionic signs.
- The space theorem at scope lines 156–169 explicitly gives the **same**
  spatial kernel's self-adjoint Hilbert–Schmidt property and
  `Z_n(beta)=sum lambda_i(beta)^n` for every even cell length. Therefore these
  are positive moments `sum |lambda_i|^n`; they give the required squaring,
  coupling, and interpolation statements without replacing the physical
  partition function by another one.
- The event measure in scope line 165 is unnormalized, of mass
  `exp(beta sum c)`. It is the Dyson-series measure, so no unreported
  partition-function normalization is introduced. Alternatively the channels
  can be normalized from norm `1/2` to one while replacing their coefficients
  `w` by `w/4`, leaving the bond unchanged. The cell rotation is unitary, and
  the ED code uses the unrotated Hamiltonian with Pauli matrices divided by
  two: diagonal entries `+/-J/4` and exchange entries `J/2`. Its energies thus
  have exactly the reported physical normalization.
- The note's squaring inequality is at
  `outputs/thermal-mpo-gap-note.tex:165–175`; coupling identities are at
  lines 218–223. Substituting separate defects `s,t` gives the asymmetric
  update printed in `ENCLOSURE_PROOF.md:182–186`, without assuming `s=t`.
- Note lines 299–313 give moment interpolation. Thus for every even
  `n<=L<=2n`, `Pt(L,beta)>=(1-s)^2(1-t)`. Note lines 119–128 and 315–329
  convert purity exceeding `1/2` into a unique ground state and the gap rate.
  Uniqueness is proved by the purity bound, not assumed from ED gaps.
- `ENCLOSURE_PROOF.md:199–220` supplies the asymmetric tail missing as a
  separately named survey theorem in this checkout. I independently checked
  its rational constants: `r_next>=r^2`, and, for `r<=1/10`,
  `r_next<5r^2`. Hence the rates decrease and their limit is at least
  `b_J-log(5)/beta_J`, precisely the core subtraction. Successive doubled
  intervals cover every even `L>=4` cells: 8,12,16,... physical spins. No claim
  is made for odd cell counts or all even physical-site counts.

## Spectral witnesses and floating-point arithmetic

The exact orbit/stabilizer enumeration is in `numerical_ed.jl:8–31`; the
verified character construction and interval block matrix are in
`verified_spectra.jl:43–85`. Characters trivial on stabilizers form an
orthonormal basis for each orbit representation; different orbits have
disjoint supports. The normalized matrix element
`sqrt(|Oj|/|Oi|) sum_y conj(chi_i(y)) H[y,rep_j]` is correct. Eighth roots of
unity and orbit-size square roots are outward enclosed, not approximated as
exact floating-point numbers. Mirroring one triangle is justified by exact
Hermiticity. Multiplicity two for charges below half filling follows from
unitary global spin flip; half filling is split into both spin-flip parities.

The assertions at `verified_spectra.jl:143–145` check every fixed-charge sum
and the full weighted dimensions: 256 at eight spins and 65536 at sixteen
spins. There are 21 and 73 nonempty solved blocks, respectively; their largest
dimensions are 14 and 1430. All these blocks use complete dense eigensolves.

I independently audited the `gamma_(4d+16)` argument, rather than taking the
stored radii as assertions of accuracy. Each real component of the residual
or Gram defect has two classical length-d dot products, at most one extra
scalar product, and at most two final additions/subtractions. This operation
budget bounds the accumulated relative errors and absolute error sources.
Fused multiply-add reduces, rather than increases, the budget. The exact
rational underflow allowance `8*2^-1022` per operation covers both gradual
underflow and FTZ/DAZ behavior. Normal-or-zero original multiplication
operands are asserted. Subnormal intermediates only enter additions or
multiplication by zero/one, so an omitted tiny value cannot undergo an
unaccounted multiplication by a large value. The overflow and `m*u<1` checks
are adequate at the actual dimensions.

The trusted arithmetic model is classical binary64 CPU OpenBLAS GEMM, using
ordinary binary64 operations/FMA, without Strassen or reduced-precision
products. This is stated, version-pinned, and backend-checked; it is not a
formal verification of OpenBLAS machine code. Julia, interval arithmetic,
transcendental rounding, LAPACK's returned finite arrays, and the runtime are
also part of the computational trust base. LAPACK accuracy is **not** trusted:
its actual residual and Gram defect are bounded afterward. The 1x1 subnormal
test is a sanity check only.

The polar-factor argument is correct: `e<1` gives an invertible Q, a unitary
U with `||Q-U||<=e` and `||Q||<2`. Since physical bond norms sum to `3L/4`,
the exact Hermitian block differs from `U Lambda U*` by at most
`r+2h+(3L/4+lmax)e`. Weyl therefore encloses all ordered eigenvalues, including
multiplicity. An isolated Ritz-value theorem is not being misused.

## Purity directions and independent executions

`verified_spectra.jl:152–161` sums every eigenvalue interval with its
multiplicity. Its binary64 shift is arbitrary and restored through enclosing
arithmetic; it need not be the exact ground energy. Lower purity bounds use
the **lower numerator** log partition and **upper denominator** log partition:
`exp(logZ_2n.lower-2logZ_n.upper)` and
`exp(logZ_2n_2beta.lower-2logZ_2n_beta.upper)`. The upper endpoint of `1-P`
is therefore an outward rational defect upper bound. `run_verified.jl`
passes these exact rationals and exact beta to the core recurrence.

The review-local script `independent_review.jl` runs the stored checker, then
regenerates all six eight-/sixteen-spin spectra for the three deltas. It
compares every stored eigenvalue center, raw residual/orthogonality component
maximum, matrix-construction error, derived norm bound, radius, symmetry
label, and multiplicity with the new calculation. Runtime fields are excluded.
Bitwise equality is required for sector labels/dimensions/multiplicities;
floating-point witness differences are counted and addressed by the rigorous
containment test below. It independently recomputes purity seeds and bounds
from the refined fresh spectral intervals. In addition,
it recomputes residuals and Gram defects with **exact complex rational matrix
multiplication for every eight-spin block**, checking squared Frobenius norms
against the independently refined bounds. This goes beyond merely reusing
stored radii.

### Independent interval containment despite floating-point differences

An initial bitwise-reproduction probe did **not** pass: 5999 checks passed
and 545 exact-equality checks failed (166.5 seconds). (Its script and log were kept on the original branch and removed
when the certificate moved to `models/`; see README.) Those checks
demanded equal approximate eigenvalues, eigenvector-dependent raw maxima,
and resulting radii/seeds; they were too strong as an acceptance condition.
The probe's printed `ACCEPTED` progress lines were unconditional progress
messages, not successful test verdicts; its final failed test summary is
authoritative. This failed reproduction attempt is not concealed.

Different binary64 results do not establish a mathematical error: classical
GEMM/eigensolves need not be bitwise reproducible, especially for degenerate
eigenvectors. The precise low-level source of each change was not traced.
Rather than assuming that these changes are harmless, the final acceptance
run independently proves containment in **every original stored eigenvalue
interval**, using a stronger forward-error bound on the fresh calculation.

Let `a=3L/4`, `h=||Htrue-H0||`'s verified bound, and use the original coarse
orthogonality bound to establish `||Q||<2`. Cauchy–Schwarz on the two
concatenated real dot-product vectors gives

```
sum_k (|Hr_ik Qr_kj| + |Hi_ik Qi_kj|) <= ||row_i(H0)|| ||col_j(Q)||
                                                   <= 2(a+h),
sum_k (|Qr_ki Qr_kj| + |Qi_ki Qi_kj|) <= 4.
```

The same inequalities hold with the real/imaginary components exchanged.
Consequently the fresh component-error allowances also admit
`gamma*(2(a+h)+qmax*lmax)+underflow` and `5*gamma+underflow`, respectively.
For each component the implementation takes the smaller of this allowance
and the original allowance: the Cauchy estimate need not improve a scalar
or very small block. It also uses the tighter exact Frobenius multiplier
`(3/2)d`, valid since `sqrt(2)<3/2`, in place of the original `2d`.
The original gamma/underflow arithmetic argument is unchanged. Reapplying
the proven polar/Weyl theorem gives a refined rigorous radius `rho_new`.
For each ordered eigenvalue the acceptance script checks, in exact rationals,

```
abs(lambda_new-lambda_stored) + rho_new <= rho_stored.
```

This directly certifies every stored spectral interval from regenerated data,
independently of its recorded eigenvector-dependent maxima. The review JSON
stores all fresh witnesses, refined radii, maximum center displacements,
positive containment margins, and production/review source hashes. Its exact
numbers are rational strings and it is written with `Results.write_result`.
Fresh refined purity defects are also checked to be no larger than the
stored rational seeds, and their computed bound is at least the stored bound.

An intermediate review-only implementation used the Cauchy bound without
taking the minimum and therefore failed containment for small blocks; that
run was interrupted (its script and log were later removed to cut size).
This was an overly large alternative error bound, not a failed enclosure of
the spectrum. The final implementation uses both rigorously valid bounds
and the improved Frobenius factor as just described. Neither exploratory
review run modified the production code or certificates.

Commands run from the repository root (originally under `certificates/`, with a
locally installed Julia 1.12.4; the paths below are the current ones):

```sh
OPENBLAS_NUM_THREADS=1 julia --project=julia models/dimerized_heisenberg/review/independent_review.jl
OPENBLAS_NUM_THREADS=1 julia --project=julia models/dimerized_heisenberg/test/verified_independent.jl
```

The four-spin analytic/exact-residual tests passed **108/108**, followed by
**2167/2167** stored witness/partition/recurrence checks. The final full-rerun
acceptance passed **4012/4012** checks in **145.5 seconds**, and successfully
wrote `acceptance.json` (about 10 MB, no longer committed;
regenerate with the first command above) with `Results.write_result`. All 63 eight-spin blocks
also passed exact rational residual/Gram checks against the refined bounds.
The final run reproduced delta=1/4 and 3/10 witnesses bitwise. Delta=1/2 had
73 blocks with changed eigenvalue arrays/radii, all independently validated
by strict interval containment; exact per-field counts are in the JSON.
Production result files were never overwritten by either review command.
A subsequent Julia `sha256_file` audit verified all six hashes in the review
JSON (three production certificates, the review script, production enclosure
source, and enclosure proof). The JSON recorded 145.782 seconds for the complete acceptance computation.

## Limitations

This accepts only the general internally reflected RP certificates, not the
original stronger real-identical-factor condition. The reference ED gaps,
bound/reference ratios, unsuccessful-grid conclusions, and correlation-length
extrapolations remain numerical diagnostics. The certified finite inputs stop
at sixteen spins; the theorem, not larger ED, supplies the stated uniform
range. The short checker validates stored witnesses and seeds but does not
reconstruct eigenvectors; this review's full reruns supply that additional
provenance check. No core or production certificate source was modified by
this review.
