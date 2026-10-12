# Verified spectral and recurrence argument

This proof applies to the **general internally reflected RP formulation** in
REPORT.md. It does not assert the impossible stronger uniform real
identical-factor identity from the original task. An independent subagent
reviewed the reflection, asymmetric recurrence, interval block construction,
and residual/polar/Weyl argument.

## Exact symmetry decomposition and block construction

At L physical spins, the Hamiltonian conserves the number q of down spins,
commutes with translation by two spins, and commutes with global spin flip.
For q<L/2 the flipped sector is unitarily equivalent and contributes a second
copy. For q=L/2 include spin-flip parity p=+/-1.

For each group orbit of representative s, a character is
`chi(t,f)=exp(2 pi i k t/(L/2))*p^f`. Keep a character precisely when it is
trivial on every stabilizer. The code tests this by the integer congruence
`2 k t + (L/2) (p==-1) f = 0 mod L`. The normalized character vector on an
orbit O is `u(g s)=chi(g)/sqrt(|O|)`. Distinct compatible characters span each
orbit's permutation representation, and disjoint orbits are orthogonal.
The fixed-charge block dimension sums and weighted total `2^L` are asserted.
These are exact representation/counting statements, not floating-point
orthogonality claims about the approximate eigensolver.

Commutation of H with the group gives the projected matrix formula

```
Hij = sqrt(|Oj|) sum_y conj(ui(y)) H[y,repj].
```

`verified_block` implements this using only the physical diagonal terms
`+/-J/4` and exchange terms `J/2`. At L=8,16 the character phases are eighth
roots of unity (L=4 is also supported). Their real/imaginary components use
only `0,+/-1,+/-sqrt(2)/2`; orbit factors are square roots of positive rational
size ratios. `IntervalArithmetic` constructs outward Float64 intervals for
all these algebraic operations and exact rational couplings.

For each upper-triangle element choose a binary64 center and compute an
outward bound a on each component of its difference from the exact element.
The exact projected matrix is Hermitian, so filling the opposite triangle by
conjugation is valid. If the block dimension is d, then
`||Htrue-H0||_2 <= ||Htrue-H0||_F <= 2 d a = h`.
The factor 2 bounds the real/imaginary Frobenius combination without needing
a rounded square root. The scalar h is stored as an exact rational.

## Rigorous residual and orthogonality bounds

A complete dense LAPACK solve produces binary64 Q and real diagonal Lambda.
No accuracy assumption about that eigensolve is made: its actual residual and
orthogonality defects are bounded afterward. Every binary64 value used in a
scalar bound is converted to its exact rational value.

Write H0=Hr+iHi and Q=Qr+iQi. The code explicitly uses real GEMMs for

```
Rr = Hr Qr - Hi Qi - Qr Lambda,
Ri = Hr Qi + Hi Qr - Qi Lambda,
Er = Qr^T Qr + Qi^T Qi - I,
Ei = Qr^T Qi - Qi^T Qr.
```

Let hmax, qmax, lmax be maximum absolute components of H0, Q, Lambda, and let
rmax, emax be maximum absolute components of the computed R and E centers.
Set `u=2^-53`, `m=4d+16`, `gamma=m u/(1-m u)`, and `tau=8*2^-1022`.
The code asserts `m u<1`, normal-or-zero original operands, finite values,
and conservative bounds excluding overflow in every product/partial sum.

The arithmetic model is standard **classical binary64 OpenBLAS real GEMM**:
binary64 multiplication/addition or fused multiply-add, with no reduced-
precision product or Strassen transformations. This is the installed CPU
OpenBLAS backend, pinned in the manifest (`OpenBLAS_jll 0.3.29+0`), with one
BLAS thread; its loaded configuration is recorded per ring. No `@fastmath`,
complex-GEMM algorithm, GPU, or mixed-precision operation is used. As with
interval elementary arithmetic, the implementation of these primitive
operations is part of the numerical trusted computing base.

For a classical dot product, chaining the elementary rounding model gives
absolute error at most `gamma_m` times the sum of absolute product terms.
There are two length-d dot products and at most one scalar product and two
final additions/subtractions in each component above. The dependency-path
and total error count are bounded conservatively by `4d+16`. Thus

```
eta_R = gamma (2 d hmax qmax + qmax lmax) + m tau/(1-m u),
eta_E = gamma (2 d qmax^2 + 1)            + m tau/(1-m u),
r = 2 d (rmax + eta_R),
e = 2 d (emax + eta_E).
```

These satisfy `||H0 Q-Q Lambda||_2<=r` and `||Q*Q-I||_2<=e`.
All displayed scalar arithmetic is implemented with `Rational{BigInt}`.

The deliberately large absolute allowance tau covers either gradual underflow
or flushing tiny intermediates to zero. Original multiplication operands are
asserted normal or zero. Tiny products/accumulators then feed only additions,
subtractions, or multiplication by zero/one; dropping subnormal inputs and
outputs costs fewer than eight smallest normal numbers per operation. There
is no unaccounted amplification of a discarded tiny multiplicative input.
The 1-by-1 subnormal BLAS check is only a sanity check; the proof does not
rely on it establishing the mode of a larger GEMM kernel.

Independent tests recompute small-block R and E in **exact complex rational
arithmetic** and compare their squared Frobenius norms with r^2 and e^2.
They also test the complete analytic four-spin spectrum and its partition
function against the generated intervals. These tests supplement, not
replace, the general forward-error argument.

## Polar-factor/Weyl enclosure theorem

Assume `e<1`, as asserted for every block. Q is invertible; its unitary polar
factor is `U=Q(Q*Q)^(-1/2)`. Since the squared singular values of Q lie in
`[1-e,1+e]`,

```
||Q-U||_2 <= e/(1+sqrt(1-e)) <= e,    ||Q||_2 < 2.
```

All physical couplings are positive for the validated parameter range and
sum to L, while `||S.S||=3/4`. Consequently `||Htrue||_2<=3L/4`. Therefore

```
||Htrue U-U Lambda||_2
 <= ||H0 Q-Q Lambda||_2 + ||(Htrue-H0)Q||_2
       + (||Htrue||_2+||Lambda||_2)||U-Q||_2
 <= r + 2h + (3L/4+lmax)e = rho.
```

Htrue and `U Lambda U*` are Hermitian. Weyl's ordered-eigenvalue inequality
then places the ith true eigenvalue in `[lambda_i-rho,lambda_i+rho]`, including
multiplicity. This is a **full-spectrum enclosure**, not an isolated
Ritz-value/residual statement that could miss levels. All centers are stored
as exact rational strings, together with rho and every scalar witness used
to recompute it.

## Partition functions and rational purity seeds

For positive rational beta, every enclosed eigenvalue is exponentiated and
summed with outward interval arithmetic, including each sector multiplicity.
An arbitrary exact binary64 energy shift stabilizes the sum:

```
log Z(beta) = -beta shift + log(sum exp[-beta(E-shift)]).
```

The shift is not assumed to be a rigorous ground energy. It is restored
exactly through enclosing arithmetic. No spectrum is truncated. The three
computed intervals for log Z_n(beta), log Z_2n(beta), log Z_2n(2beta) give

```
Ps >= exp(logZ_2n.lower - 2 logZ_n.upper),
Pt >= exp(logZ_2n_2beta.lower - 2 logZ_2n_beta.upper).
```

The code uses interval expressions for these same formulas, then takes the
upper endpoint of `1-P`. That endpoint is an exact binary64 rational; it is
saved as a rational string and passed with exact beta to
`certified_gap_bound`. The checker reconstructs all log-partition intervals
from stored spectral centers/radii and recomputes the seeds and bound.
It does not rerun the expensive eigenvector solves; reproducing those is the
purpose of `run_verified.jl`.

## General reflected RP and the asymmetric recurrence

The cell reflection proof in REPORT.md establishes a real symmetric bounded
Hilbert--Schmidt spatial kernel, hence the even moment representation
`Z_n(beta)=sum |lambda_i(beta)|^n`. Internal reflection is an algebra
automorphism on operators, so it preserves chronological product order; no
unproved time-order reversal or fermion sign is involved.

The survey theorem numbers mentioned in `certified.jl` are not separately
available as named proofs in this checkout. The following derivation instead
uses explicit proved identities of `outputs/thermal-mpo-gap-note.tex`, extended
to this moment form by Appendix E of the scope note.

Let `s,t` upper-bound the two defects and put

```
z(s,t)=1-(1-s)^2(1-t),   f(z)=z^2/[2(1-z)^2].
```

Equation `couple1`, using thermal purity <=1, gives
`Ps(n,2beta)>=(1-s)^2(1-t)=1-z(s,t)`. Equation `square` gives
`s'=f(z(s,t))` at `(2n,2beta)`. Equation `couple2`, using spatial purity <=1,
gives `Pt(4n,beta)>=(1-t)^2(1-s')`; another `square` gives
`t'=f(z(t,s'))`. This is precisely the implemented unequal-defect map.

For every even `n<=L<=2n`, equations `moments` and `interpolate` give

```
Pt(L,beta) >= [(1-t)(1-s)^2]^(L/(2n)) >= 1-z(s,t).
```

If `z<1/2`, equations `weight` and `gapratio` imply a unique ground state and
`Delta_L >= -log(r)/beta`, where `r=z/(1-z)`. Every even L>=n0 belongs to a
successive doubled interval. No separate uniqueness assumption is borrowed
from a numerical reference gap.

For completeness, the rigorous tail used by the core follows as well.
Write `beta_j=2^j beta0`, `b_j=-log(r_j)/beta_j`. Since `s'=r^2/2` and t'>=0,

```
r_next >= (1-s')^-2-1 >= r^2,
```

so b_j is nonincreasing. For r<=1/10, t<=r and
`z(t,s')<=2r+r^2/2<=(41/20)r`. It follows exactly that

```
t' <= (84050/25281) r^2,
z_next <= (109331/25281) r^2,
r_next <= (10933100/2418769) r^2 < 5 r^2.
```

This preserves r<=1/10. Iterating the upper inequality and summing the
geometric tail gives

```
b_infinity >= b_J - log(5)/beta_J.
```

Thus the tail subtraction and the uniform bound returned by
`certified_gap_bound` are justified under the general reflected RP hypothesis.
Its outward interval orbit also checks the necessary domain/admissibility
conditions. The original stronger real-identical-factor condition remains
false; it is not used in this general-RP certificate.
