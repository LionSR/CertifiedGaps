# Method: verified spectral, recurrence and exclusion arguments

This proof applies to the **general internally reflected RP formulation** of README.md
(`-h = sum C tensor conj(C)` in the cell basis below, complex channels allowed). It does not
assert the stronger uniform real identical-factor identity, which is false for this bond.
All theorem inputs are exact rationals or outward intervals; LAPACK and KrylovKit outputs are
never accepted as inequalities. Entry points: `certify.jl`, `check.jl`, `src/Certification.jl`.

## Exact reflection positivity

Choose strong pairs as cells. Set `s=1+delta`, `w=1-delta`, `A_a=S_a tensor I`,
`B_a=I tensor S_a`. Then the on-site term is `D=s sum A_a B_a` and the bond is
`h=w sum B_a tensor A_a`. Let

```
J = [0 1; -1 0], R = I tensor J,
M = [1  1  0   0;
     0  0  1   i;
     0  0  1  -i;
     1 -1  0   0],    U = R M/sqrt(2).
```

`M/sqrt(2)` consists of three real triplet combinations and the singlet multiplied by i;
`M'M=2I` and `M M^T/2=F`, the internal cell swap. In the rotated product frame

```
Drot = s(-A_x B_x + A_y B_y - A_z B_z),
-hrot = w(B_x tensor A_x - B_y tensor A_y + B_z tensor A_z).
```

The antiunitary `F K` exchanges the two cell ends, fixes `Drot`, and maps the real channels
`(B_x, iB_y, B_z)` to `(A_x, iA_y, A_z)`; this is the antiunitary bond form of Appendix E of
`notes/thermal-mpo-gap-scope.tex`. Changing to `M/sqrt(2)` turns the antiunitary into ordinary
conjugation. Explicitly, with `C_a = sqrt(w) K_a/2`,

```
K_x = [0 0 1 0; 0 0 0 i; 1 0 0 0; 0 -i 0 0],
K_y = [0 0 0 i; 0 0 1 0; 0 -1 0 0; i 0 0 0],
K_z = [0 1 0 0; 1 0 0 0; 0 0 0 -i; 0 0 i 0],
```

exactly `-(U tensor U)' h (U tensor U) = sum C_a tensor conj(C_a)` and
`U' D U = s diag(-3,1,1,1)/4`, which is real. The two-spin cell has one internal bond, so the
palindromic condition is automatic. `src/rp.jl` checks these identities in
`Complex{Rational{BigInt}}`, cancelling the normalisation denominators symbolically. Each seed
certificate stores the model and the factors as a core `RingModels` record (rational weight `w`,
matrices `K_a/2`); the checker rebuilds the model from `delta` and verifies the stored factors
with `RingModels.verify_rp`.

## Exact symmetry decomposition and block construction

At L physical spins, the Hamiltonian conserves the number q of down spins, commutes with
translation by two spins, and commutes with global spin flip. For q<L/2 the flipped sector is
unitarily equivalent and contributes a second copy. For q=L/2 include spin-flip parity p=+/-1.

For each group orbit of representative s, a character is
`chi(t,f)=exp(2 pi i k t/(L/2))*p^f`. Keep a character precisely when it is trivial on every
stabilizer; the code tests the integer congruence `2 k t + (L/2) (p==-1) f = 0 mod L`. The
normalized character vector on an orbit O is `u(g s)=chi(g)/sqrt(|O|)`. Distinct compatible
characters span each orbit's permutation representation, and disjoint orbits are orthogonal.
Each charge's block dimensions must sum to `binomial(L,q)` and the weighted total to `2^L`;
these exact counting checks also catch a missing momentum or parity block.

Commutation of H with the group gives the projected matrix formula

```
Hij = sqrt(|Oj|) sum_y conj(ui(y)) H[y,repj].
```

`verified_block` implements this using only the physical diagonal terms `+/-J/4` and exchange
terms `J/2`. At L=4,8,16 the character phases are eighth roots of unity, with components
`0,+/-1,+/-sqrt(2)/2`; orbit factors are square roots of positive rational size ratios.
`IntervalArithmetic` encloses all these algebraic operations and the exact rational couplings.
For each upper-triangle element choose a binary64 centre and an outward bound a on each
component of its error; exact Hermiticity justifies filling the lower triangle by conjugation.
Then `||Htrue-H0||_2 <= ||Htrue-H0||_F <= 2 d a = h`, stored as an exact rational.

**Conjugate momenta.** The physical H is real, so character conjugation gives
`H_(n-k,p) = conj(H_(k,p))` in the same orbit ordering (compatible stabilizers agree). For
k>n/2 the spectrum and radius of `(n-k,p)` are reused; this is separate from the spin-flip
multiplicity. Every sector still appears in the certificate and in the dimension accounting;
the fast checker requires a reused block to repeat its partner exactly, and the tests compare
independently constructed conjugate interval blocks for equality.

## Rigorous residual and orthogonality bounds

A complete dense LAPACK solve produces binary64 Q and real diagonal Lambda. No accuracy
assumption about that eigensolve is made: its residual and orthogonality defects are bounded
afterwards. Every binary64 value used in a scalar bound is converted to its exact rational.

Write H0=Hr+iHi and Q=Qr+iQi. The code explicitly uses real GEMMs for

```
Rr = Hr Qr - Hi Qi - Qr Lambda,
Ri = Hr Qi + Hi Qr - Qi Lambda,
Er = Qr^T Qr + Qi^T Qi - I,
Ei = Qr^T Qi - Qi^T Qr.
```

Let hmax, qmax, lmax be the maximum absolute components of H0, Q, Lambda, and rmax, emax those
of the computed R and E. Set `u=2^-53`, `m=4d+16`, `gamma=m u/(1-m u)`, `tau=8*2^-1022` and
`eta = m tau/(1-m u)`. The code asserts `m u<1`, normal-or-zero original operands, finite
values, and conservative bounds excluding overflow in every product and partial sum.

The arithmetic model is **classical binary64 OpenBLAS real GEMM**: binary64 multiplication and
addition or fused multiply-add, with no reduced-precision product or Strassen transformation.
This is the installed CPU OpenBLAS backend, pinned in the manifest (`OpenBLAS_jll 0.3.29+0`),
with one BLAS thread; its loaded configuration is recorded per ring. No `@fastmath`,
complex-GEMM algorithm, GPU or mixed-precision operation is used. As with interval elementary
arithmetic, these primitive operations are part of the numerical trusted computing base. For a
classical dot product, chaining the elementary rounding model gives absolute error at most
`gamma` times the sum of absolute product terms, plus `eta`; the dependency paths contain fewer
than m elementary operations. The allowance `tau` covers gradual underflow or flushing tiny
intermediates to zero: original multiplication operands are asserted normal or zero, so tiny
products and accumulators feed only additions, subtractions, or multiplication by zero or one.

**Two passes.** First, using only component maxima,

```
e_coarse = 2d [emax + gamma(2 d qmax^2 + 1) + eta] < 1
```

bounds `||Q*Q-I||_2` (the factor 2d bounds the 2-norm of a complex d-by-d matrix by its
component maximum). This proves every column of Q has Euclidean norm below sqrt(2). By
Cauchy-Schwarz on the concatenated real vectors, each absolute Gram dot-product sum is then at
most 2 (plus 1 for the identity), and each absolute residual dot-product sum at most
`2||H0|| <= 2(3L/4+h)` (the couplings are positive and sum to L, `||S.S||=3/4`). The sharper,
non-circular second pass is

```
r = 2d [rmax + gamma(2(3L/4+h) + qmax lmax) + eta],
e = 2d [emax + 3 gamma + eta] < 1,
```

with `||H0 Q-Q Lambda||_2 <= r` and `||Q*Q-I||_2 <= e`. All of this scalar arithmetic is exact
`Rational{BigInt}`. Tests recompute R and E in exact complex rational arithmetic for every
4-spin block (six deltas) and every 8-spin block (three deltas) and compare their squared
Frobenius norms with r^2 and e^2; they supplement, not replace, the forward-error argument.

## Polar-factor/Weyl enclosure theorem

Assume `e<1`. Q is invertible; its unitary polar factor is `U=Q(Q*Q)^(-1/2)`. Since the squared
singular values of Q lie in `[1-e,1+e]`,

```
||Q-U||_2 <= e/(1+sqrt(1-e)) <= e,    ||Q||_2 < 2.
```

With `||Htrue||_2 <= 3L/4`,

```
||Htrue U-U Lambda||_2
 <= ||H0 Q-Q Lambda||_2 + ||(Htrue-H0)Q||_2 + (||Htrue||_2+||Lambda||_2)||U-Q||_2
 <= r + 2h + (3L/4+lmax)e = rho.
```

Htrue and `U Lambda U*` are Hermitian. Weyl's ordered-eigenvalue inequality places the ith true
eigenvalue in `[lambda_i-rho,lambda_i+rho]`, including multiplicity. This is a **full-spectrum
enclosure**, not an isolated Ritz-value statement that could miss levels. Tests compare the
entire four-spin spectrum with its closed form: singlets `-1 +/- sqrt(1+3delta^2)`, triplets
`-1, -delta, +delta` (three each) and the quintet `1`.

## Partition functions and rational purity seeds

For exact positive beta, every enclosed eigenvalue is exponentiated and summed with outward
interval arithmetic, including each sector multiplicity. An arbitrary exact binary64 energy
shift stabilises the sum: `log Z(beta) = -beta shift + log(sum exp[-beta(E-shift)])`. The shift
is not assumed to be a ground energy, and no level is discarded. The three intervals for
`log Z_n(beta)`, `log Z_2n(beta)`, `log Z_2n(2beta)` give

```
Ps >= exp(logZ_2n.lower - 2 logZ_n.upper),
Pt >= exp(logZ_2n_2beta.lower - 2 logZ_2n_beta.upper).
```

The upper endpoints of `1-P` are rounded up to the decimal grid `10^-6` plus one spare unit and
stored as exact rationals `eps_s`, `eps_t`; the spare unit lets an independent replay on another
BLAS prove the same stored seeds. `certified_gap_bound(eps_s, eps_t, beta)` gives the bound,
whose exact binary value is stored and which is displayed rounded down to `10^-6`. For n0=4
cells the physical range is **N divisible by 4 and N >= 8**, not every even site count.

The fast checker re-derives each block's witness and radius from its stored component maxima,
and recomputes the log-partition intervals, seeds and bound from the stored centres. `check.jl
--full` trusts neither: it rebuilds every block, eigenvector, residual and construction error,
requires the fresh intervals to meet the stored ones, and proves the stored seeds, the 16-spin
gap interval and every exclusion box again from the fresh spectra only.

## General reflected RP and the asymmetric recurrence

The cell reflection establishes a real symmetric bounded Hilbert--Schmidt spatial kernel, hence
the even moment representation `Z_n(beta)=sum |lambda_i(beta)|^n`. Internal reflection is an
algebra automorphism on operators, so it preserves chronological product order; no unproved
time-order reversal or fermion sign is involved. The following derivation uses explicit proved
identities of `notes/thermal-mpo-gap-note.tex`, extended to this moment form by Appendix E of
the scope note.

Let `s,t` upper-bound the two defects and put

```
z(s,t)=1-(1-s)^2(1-t),   f(z)=z^2/[2(1-z)^2].
```

Equation `couple1`, using thermal purity <=1, gives `Ps(n,2beta)>=(1-s)^2(1-t)=1-z(s,t)`.
Equation `square` gives `s'=f(z(s,t))` at `(2n,2beta)`. Equation `couple2`, using spatial
purity <=1, gives `Pt(4n,beta)>=(1-t)^2(1-s')`; another `square` gives `t'=f(z(t,s'))`. This is
precisely the unequal-defect map of the core.

For every even `n<=L<=2n`, equations `moments` and `interpolate` give

```
Pt(L,beta) >= [(1-t)(1-s)^2]^(L/(2n)) >= 1-z(s,t).
```

If `z<1/2`, equations `weight` and `gapratio` imply a unique ground state and
`Delta_L >= -log(r)/beta`, where `r=z/(1-z)`. Every even L>=n0 belongs to a successive doubled
interval. No uniqueness assumption is borrowed from a numerical reference gap.

The rigorous tail used by the core: write `beta_j=2^j beta0`, `b_j=-log(r_j)/beta_j`. Since
`s'=r^2/2` and t'>=0, `r_next >= (1-s')^-2-1 >= r^2`, so b_j is nonincreasing. For r<=1/10,
t<=r and `z(t,s')<=2r+r^2/2<=(41/20)r`. It follows exactly that

```
t' <= (84050/25281) r^2,
z_next <= (109331/25281) r^2,
r_next <= (10933100/2418769) r^2 < 5 r^2.
```

This preserves r<=1/10. Iterating and summing the geometric tail gives
`b_infinity >= b_J - log(5)/beta_J`. Thus the tail subtraction and the uniform bound returned
by `certified_gap_bound` are justified under the general reflected RP hypothesis; its outward
interval orbit also checks the domain and admissibility conditions.

## Rigorous negative results for every beta > 0

A beta-grid failure alone would not establish the grid minimum. `negative_cover` proves that
the recurrence above fails for every beta at n0=2 cells and, where no seed exists, n0=4:

1. Thermal purity is nondecreasing in beta, since its logarithmic derivative is
   `2(<H>_beta - <H>_2beta) >= 0`. A verified `P_tau(2n,1) <= 1/2` (thermal defect at least
   1/2, hence z >= 1/2) excludes `0 < beta <= 1`.
2. If `dE = E0(2n)_lower - 2 E0(n)_upper > 0`, positivity and dimension counting give
   `P_s(n,beta) <= 2^(4n) exp(-beta dE)`. An exact dyadic B making this `<= 1/2` excludes
   every `beta >= B`.
3. `[1,B]` is subdivided. On each closed box, interval evaluation of the shifted partition sums
   bounds both true defects from below. These optimistic lower bounds are propagated through
   `z`, `s' = f(z)` and the analogous t update, rounding **down** to `10^-15` after each step;
   the map is monotone, so lower bounds stay lower bounds. Once `z >= 1/2` the core recurrence
   fails (it requires `z < 1/2` at every step); if the second coupling `zt >= 1/2` then
   `t' >= f(1/2) = 1/2` and the next z also reaches the threshold.

The stored boxes are contiguous and exhaust `[1,B]`; the checker re-evaluates every box and
both tails. These negatives exclude **this recurrence at these seeds** (n0 = 2, 4, the only
even seeds whose doubled ring has at most 16 spins). They do not exclude larger rings, other
seed theorems (for example the survey's Theorem 2.24(v) replacement), other methods, or a
physical gap.

## A scalable route beyond ED: dimer perturbation

At w=0, `Z_n=[exp(3 beta s/4)+3 exp(-beta s/4)]^n`. For w>0 write `H=H_dimer+V` with
`||V_n|| <= v_n = 3nw/4`; Weyl gives `|log Z_n - log Z_n^dimer| <= beta v_n`. Hence, for any n,

```
P_s(n,beta) >= exp(-4 beta v_n),
P_tau(2n,beta) >= [(1+3x^2)/(1+3x)^2]^(2n) exp(-8 beta v_n),   x = exp(-beta s).
```

These bounds cost O(log n) operations, independent of the Hilbert dimension. The tests close
an actual seed at delta=999999/1000000, n0=128, beta=8 (512 spins in the largest ring). At the
requested deltas even n0=4, beta=6 gives a spatial-defect upper bound above 1/2, so the
estimate certifies nothing there, and increasing n worsens the extensive error. A linked
cluster expansion that keeps the bulk free-energy cancellation would be needed.

## Floating-point diagnostics

`scripts/scan20.jl` enumerates every Sz, momentum and spin-flip block of the 20-spin ring.
Large blocks keep 40 converged Krylov Ritz values. The omitted-level estimate assumes these are
all the lowest levels with correct multiplicities, which is not verified. All 20-spin gaps,
thermal-tail estimates, odd n0=5 candidates and energy-decay fits are numerical and never enter
the positive or negative certificates.
