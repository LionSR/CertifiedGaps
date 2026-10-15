# Heisenberg ladder: fresh certificate and threshold audit

The current machine-checkable result is `results/certified_suite.json`, produced
by `certify_suite.jl` and checked from fresh complete spectra by `check.jl`.
The original J=4 certificate and numerical scans remain under `results/` as
historical evidence. The README reports the current bounds and run times.

## Hamiltonian and length convention

H_L = sum_{j=1}^L sum_{a=1}^2 S_{a,j}.S_{a,j+1}
      + J sum_{j=1}^L S_{1,j}.S_{2,j}, with periodic j, spin 1/2 and leg coupling 1.
One cell is a rung, d=4; there are N=2L spin sites. At L=2 the leg bonds are
counted twice, as in the repository's ring convention. At even L the rotated
Hamiltonian is translation invariant with one rung per cell.

The criterion requires even n0. With rings no larger than ten rungs, direct
seeds have n0=2 or 4. Odd rings in the historical scans are diagnostics, and
their n0=3,5 entries are excluded from every claim. A certificate at n0=4
covers all N=8,12,16,... spin sites, not every integer N>=8.

## Exact reflection positivity and equivalence

Use the exact rational matrices

    X = [0 1; 1 0]/2, A = i Sy = [0 1; -1 0]/2,
    Z = [1 0; 0 -1]/2, u = 2A = exp(i pi Sy).

The original two-spin bond h=X tensor X - A tensor A + Z tensor Z becomes
h'=-X tensor X - A tensor A - Z tensor Z after rotation of exactly one endpoint.
For even L the given checkerboard sublattice rotates one endpoint of every
rung and leg bond. On alternating rungs the rotations are u tensor I and
I tensor u. Both give the same real symmetric on-rung term J h'. Both bond
orientations give

    -h_cell = sum_C C tensor C,
    C in {X tensor I, A tensor I, Z tensor I,
          I tensor X, I tensor A, I tensor Z}.

`rp.jl` verifies the original-to-rotated identities over Rational{BigInt};
`rp_model` also stores the explicit matrices and uses the shared core's exact
RP checker. This justifies computing spectra in the unrotated SU(2) basis.
No floating-point RP screen is accepted as a proof.

## Complete verified spectra

The original spin basis is resolved by down-spin number m, rung momentum k,
leg exchange and (at m=L) spin inversion. Orbits and stabilizers are exact
permutations. A character is included exactly when it is trivial on the
stabilizer. Opposite magnetizations and conjugate momenta have exact spectral
multiplicities. The dimension sum must be 4^L.

The projected entries are rational coefficients times enclosed roots of unity
and square roots of exact orbit-size ratios. An interval row-sum bound herror
controls the distance from the exact Hermitian block to its Float64 center H.
LAPACK supplies tentative Q,Lambda. The implementation bounds

    r >= ||H Q - Q Lambda||,
    eta >= ||Q* Q - I|| < 1,
    d = eta/(1+sqrt(1-eta)).

Frobenius norms use interval accumulation. Complex dot products receive the
explicit error allowance gamma_(8n+32), with unit roundoff 2^-53, and a 2^-1000
underflow allowance. The latter exceeds all possible underflow accumulation
for the block dimensions used here. The polar factor U of Q satisfies
||Q-U|| <= d, so Weyl gives a simultaneous sorted-eigenvalue radius

    rho = r + (||H|| + max |Lambda|) d + herror.

Every eigenvalue and multiplicity enters the outward-rounded exponential sums.
Thus the certificate encloses log Z_n(beta), log Z_2n(beta), log Z_2n(2 beta)
and both defects. Upper defects are rounded upward to the exact 10^-8 grid.
The shared `certified_gap_bound` checks the recurrence and its infinite tail;
the reported gap is rounded downward to an exact multiple of 10^-9.
The machine assumptions are IEEE binary64 round-to-nearest for OpenBLAS with
one BLAS thread and the repository's abstract criterion.

## Threshold: what is certified

The suite re-evaluates the entire grid n0 in {2,4}, beta=k/4 for 1<=k<=64 at
all six requested couplings and at the two rational endpoints proposed by the
old numerical threshold search. Successes use upper defects. Failures use
lower defects and a separate exact-rational witness; a failed upper-bound
attempt alone is recorded as unresolved.

For lower defects s,t, put z=1-(1-t)(1-s)^2 and f(z)=z^2/[2(1-z)^2]. If z>=1/2,
the next spatial and thermal components are both >=1/2. If the thermal
intermediate argument reaches 1/2 first, the next thermal component is >=1/2
and the following spatial argument reaches 1/2. The order-preserving
recurrence therefore cannot converge to zero. Intermediate lower bounds are
rounded downward to 10^-24, which preserves this rejection proof.

A passing coupling next to a grid-rejected coupling locates a **local grid
transition**, without assuming monotonicity in J. It is not a global infimum
over continuous beta, larger seeds, other methods, or all couplings. It is
not the physical gap-closing coupling. The code never turns the old numerical
bisection or a finite failed scan into any of those stronger statements.

## Continuous-temperature exclusion at J=2 and J=3/2

For each eligible n, the complete enclosures give a_L<=E0(L)<=b_L. Set
F_L(beta)=log Tr exp[-beta(H_L-a_L)]. Since H_L-a_L is positive semidefinite,
F_L is decreasing. For beta in [l,u],

    log P_s(n,beta) <= F_2n(l) - 2 F_n(u) - beta (a_2n-2a_n).

The right side is evaluated as an interval, so no spatial-purity monotonicity
is assumed. Thermal purity is increasing: its log derivative is
2(<H>_beta-<H>_2beta)>=0, because d<H>_beta/d beta=-Var_beta(H)<=0.
Its value at u therefore bounds P_tau(2n,beta) from above. These give lower
seed defects on a whole temperature interval; exact recurrence rejection
certifies that entire interval.

For the infinite-beta tail, delta=a_2n-2b_n>0 and, for beta>=B,

    log P_s(n,beta) <= F_2n(B) - beta delta <= F_2n(B) - B delta.

The checker verifies positive delta, rejection at B, and contiguous coverage
from beta=0 to B by stored interval witnesses. A successful exclusion rules
out these two seeds for every beta>0, but says nothing about larger seeds or
the physical spectral gap.

## Routes beyond ED

`IDEAS.md` compares closed forms, convergent expansions, transfer contraction,
perturbation, smaller-seed identities, symmetry counting and variational energy
sandwiches. `PERTURBATION.md` proves that even a centered independent-rung norm
sandwich fails for every seed length and beta when 0<J<=4. That arbitrary-size
route has a clean negative result; no unproved cluster remainder is claimed.

The same norm argument is useful locally in coupling. For |J-J*|<=delta,
centering the rung perturbation gives ||W_L+LI/4||<=L/2 and the rigorous seeds

    s_delta = 1-(1-s) exp(-2 beta n delta),
    t_delta = 1-(1-t) exp(-4 beta n delta).

Whenever their recurrence closes, the suite stores a uniform bound on the
entire coupling interval. This requires no monotonicity in J.

## Historical scan audit and reproducibility

`replay_scans.jl` re-derives every symmetry block dimension through L=10 and
compares the exact label/multiplicity maps with each of the six stored dense
scans. It reproduces all 256 grid rows per coupling from the retained log Z
grids and records SHA-256 hashes of every source artifact. It does not rerun
those long dense diagonalizations; their ten-rung numbers remain floating point.
The fresh certificate uses no historical numerical spectrum or seed as a
trusted input. `results/scan_replay.json` records this distinction explicitly.

The model has its own Julia 1.12 project, a path dependency on `../../julia`,
and a pinned manifest. Every executable source and the model environment are
hashed in each new certificate. `check.jl` checks the hashes, exact RP,
fresh complete spectra, seed domination, recurrence, all negative witnesses,
coupling continuation and the claimed length range. Tests additionally compare
small sectors with an independently built full graph Hamiltonian and reject
corrupted RP, missing temperature coverage and an invalid tail witness.
No core or other model files were changed.

For planning only, the retained floating-point extrapolation in
`results/extrapolation.json` suggests n0=6 rungs at J=2 and n0=8 at J=3/2,
requiring 12- and 16-rung traces (24 and 32 spins). It fits an exponential to
finite-size ground-energy mismatches and uses a dilute-triplet thermal model.
These estimates have no controlled remainder and are not seed enclosures or
certified correlation lengths; no proof in this report relies on them.
