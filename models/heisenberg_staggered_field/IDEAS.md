# Candidate routes

The rotated frame has a one-site unit cell, despite the two-site period of the physical
field. Only even lengths are equivalent by the sublattice rotation. Thus n0 must be even
and counts individual spins. Translation by two sites and total Sz commute in the physical
frame and give momentum/charge blocks.

1. Dense symmetry blocks: dimension approximately binomial(L,L/2)/(L/2), cubic eigensolve
   cost per block and quadratic storage. This is practical for the stated h=1 seed (8,5.25),
   which requires L=16, but L=20/24 become substantially more expensive. Chosen route.
   Verify eigenvalues using residual and orthogonality bounds, including floating-point
   matrix-product error and errors in the momentum-block entries.
2. Large-field perturbation: could bound partition functions independently of ED size,
   but h<=1 is not a small perturbation parameter for exchange strength one.
3. Cluster expansion: potentially gives large rings but no explicit convergence/remainder
   estimate is supplied for these fields; developing one is beyond the finite ED search.
4. MPS/DMRG or partial Krylov spectra: useful for numerical exploration only; cannot
   replace complete dense sector spectra on the certificate path.
5. Smaller-seed variant: not selected because the requested criterion explicitly uses
   the thermal purity of the doubled ring.

First verify RP exactly and the environment. Search numerical spectra on a beta grid,
then rerun the best feasible seed with enclosures. If computational resources limit a
search, distinguish that from an actual absence of seeds or a physical gap threshold.

## Outcome

Route 1 certified all five fields at n0 = 8 from complete 8- and 16-site spectra. The best
20-site candidates (n0 = 10) were found with partial Krylov spectra and remain uncertified:
the largest 20-site block (dimension 18504) is beyond the dense verifier's memory budget. The
shared-core cross-check (`scripts/core_certificate.jl`) reproduces the h = 1 seed with a
different block structure and the exact rational chain.
