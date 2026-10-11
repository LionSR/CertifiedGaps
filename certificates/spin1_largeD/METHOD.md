# Model, sectors and enclosure proof

## Exact model and reflection positivity

In the Cartesian basis take the integer matrices

```
Gx = [0 0 0; 0 0 1; 0 -1 0]
Gy = [0 0 -1; 0 0 0; 1 0 0]
Gz = [0 1 0; -1 0 0; 0 0 0].
```

Set Sa = -i Ga. The checker verifies [Sx,Sy]=i Sz and sum Sa^2=2I exactly.
The independently constructed spin-one bond is h=P-|Omega><Omega|, where P swaps
factors and Omega=|xx>+|yy>+|zz>. Integer arithmetic verifies

    -h = Gx⊗Gx + Gy⊗Gy + Gz⊗Gz.

Thus every coefficient is exactly 1 and every G is real. The onsite term is
D diag(1,1,0)=D Sz², real symmetric, for every rational D in the scan.
This establishes the local RP condition in the supplied scope note. It does not
by itself establish the seed hypothesis.

In the normalized Sz basis m=-1,0,1, each allowed nearest-neighbour transition
(mj,mj+1) -> (mj+1,mj+1-1) has amplitude 1: the two sqrt(2) ladder factors
multiply and the bond contributes 1/2. The diagonal is
sum mj*mj+1 + D sum mj². These are exact rational numbers; L=2 has two directed
bonds, as required by the periodic Hamiltonian convention in the core.

## Symmetry accounting

All product states are partitioned by total Sz. Charge q and -q have identical
spectra because m -> -m commutes with H. In q=0, use translation and spin
inversion together, with character exp(2πikt/L) p^f, p=±1. In nonzero charge,
use translation alone. Orbits contribute only when the character is trivial
on their stabilizer. Each projected orbit vector is normalized by sqrt(orbit
size). A transition from representative a to T^u F^f b contributes

    h * sqrt(size(a)/size(b)) * exp(-2πiku/L) * p^f.

Distinct projected orbits are orthonormal. The group commutes exactly with H,
so the blocks are invariant. Complex conjugation pairs k and L-k. Each computed
block therefore has charge multiplicity 1 or 2 and momentum multiplicity 1 or 2.
The explicit assertion sum(multiplicity*dimension)==3^L prevents missing sectors.
The certification path uses dense Hermitian eigensolvers for every block; no
partial spectrum enters a certified partition function. Small tests independently
assemble the full product-space H and compare all eigenvalues through L=6.

## Matrix and eigenvalue enclosures

Every rational transition, trigonometric phase and square root in the block
construction is evaluated in outward-rounded IntervalArithmetic. The diagonal
and the independent upper triangle determine a Hermitian matrix A. Interval
entrywise deviations from A give a maximum absolute row-sum bound c on
||H_block-A||₂ (Hermiticity also bounds the column sum). Physical operator norms
give ||H_block||₂ <= L(2+|D|), hence ||A||₂ <= M=L(2+|D|)+c.

Let A Q ≈ Q Lambda be the full dense eigensystem, with real diagonal Lambda.
We actually compute R_fl=fl(AQ-QLambda) and E_fl=fl(Q*Q-I), where * is the
conjugate transpose. All block dimensions and eigenvalue multiplicities are
retained. For a block of dimension m, set u=2^-53 and

    gamma = (8m+32)u / (1-(8m+32)u).

The standard componentwise floating-point dot-product error model implies

    ||AQ-QLambda||₂ <= r
      = ||R_fl||F + gamma (sqrt(m) M ||Q||F + ||Q||F ||Lambda||₂) + tiny,
    ||Q*Q-I||₂ <= e
      = ||E_fl||F + gamma (||Q||F² + sqrt(m)) + tiny.

The 8m+32 factor conservatively covers real operations in complex products,
accumulation and subtraction. Every norm upper bound and scalar expression is
outward rounded. A floating sum of squared magnitudes of N entries is increased
using gamma_(4N+16), then square-rooted outward. The additive 10^-290 inside a
squared norm and 10^-280 in residual bounds cover subnormal rounding effects
by many orders of magnitude at these dimensions. The implementation assumes
IEEE binary64 arithmetic and the standard BLAS dot-product error model; it
checks e<1 and runs with one OpenBLAS thread. It does not assume that LAPACK's
returned eigenvalues are exact or invoke its undocumented backward error.

Here is the self-contained polar/Weyl theorem used. Because e<1, Q is invertible.
Its polar factor U is unitary, and

    ||Q-U||₂ <= d = e/(1+sqrt(1-e)).

Therefore

    ||AU-U Lambda||₂ <= r + (M+||Lambda||₂)d.

Multiplication on the right by U* gives the same bound on the difference of the
Hermitian matrices A and U Lambda U*. Weyl's eigenvalue perturbation inequality
then bounds each sorted eigenvalue, with multiplicity, by the corresponding
sorted diagonal entry of Lambda plus/minus

    eta = c + r + (M+||Lambda||₂)d.

Weyl itself follows directly from the min-max characterization: adding a
Hermitian perturbation of norm eta shifts every Rayleigh quotient by at most
eta, hence shifts every ordered min-max value by at most eta. This argument
also proves the enclosure in degenerate eigenspaces.

## Partition functions and seeds

For the complete spectrum, sum exp(-beta*[lambda-eta,lambda+eta]) with outward
rounding after an exact binary energy shift, then take log and undo the shift.
Each partition-function interval contains the true log Z. Use

    eps_s = upper(1-exp(log Z_2n(beta)-2 log Z_n(beta))),
    eps_t = upper(1-exp(log Z_2n(2 beta)-2 log Z_2n(beta))).

The binary upper endpoints are converted exactly to Rational{BigInt}, stored as
strings, and passed with rational beta to the unmodified certified_gap_bound.
The seed assembler also applies the common maximum eigenvalue radius to the
first two sorted levels of each seed ring and asserts a strictly positive
finite-ring gap, verifying ground-state uniqueness on those two rings.
The lower result is the uniform gap certificate under the repository's supplied
two-purity theorem, for every even L>=n0 one-site cells, hence sites. Only even
n0 are eligible; odd-n scan values are diagnostics and are never promoted to
certificates.

The fast checker reruns exact RP, verifies source/data fingerprints, recomputes
partition intervals from stored eigenvalue enclosures, and recomputes the
recurrence from stored rational seeds. Rebuilding and verifying all dense
eigensystems requires certify.jl and is intentionally outside the fast checker.
