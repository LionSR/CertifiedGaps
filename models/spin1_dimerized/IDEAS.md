# Candidate routes

The authorized numerical route is full dense ED in Sz, cell momentum and spin-flip
sectors, followed by verified residual enclosures only at the best eligible seed.
At twelve spins the largest unreduced Sz sector is 73789; cell momentum and spin flip
reduce it, but full dense spectra remain expensive (cubic work, quadratic memory).
The six-cell cap permits n=2,3 only; the repository theorem requires even n, hence
only n=2 is an eligible certificate seed. n=3 is a numerical diagnostic, not a
uniform-gap certificate under the stated theorem. n=4 would need sixteen spins.

A strong-dimer perturbation expansion could bypass ED size limits, but implementing
and checking an explicit convergent remainder is a separate substantial project.
The exact decoupled dimer point is a useful validation case, not one of the requested
parameters. A transfer-kernel argument using swap reflection inside the cell is
needed: the intercell interaction is not literally a sum of identical cell operators.

Chosen route: first verify site RP and the cell-reflection kernel identity; implement
symmetry-block dense ED with exact integer spin-1 Hamiltonian entries times rational
couplings. Use a rigorous residual/orthogonality norm bound for the best even seed,
and outward interval partition sums. Report any computational limits honestly.

## Refinement after inspecting the event expansion

The initial even-cell restriction is unnecessarily strong for this model. Write the
one-site real event kernel as F(a,b). It is symmetric because the insertions from
the two incident bonds are the same real matrices. With positive weak/strong event
measures the blocked kernel is X(a,b)=integral F(a,u)F(b,u) dmu_weak(u)=KK*.
Thus X is positive semidefinite, not merely self-adjoint, and the proof of the
criterion applies also to odd n=3 cells: positivity of X^n and the moment
interpolation are the only uses of even n in the supplied proof. We retain the
requested conclusion for even cell lengths. The six-cell cap still excludes n>=4.

For rigorous twelve-spin blocks, use real dihedral projectors (cell momentum plus
reflection) and spin flip. For cell periods 2,3,4,6, twice every cosine in the real
character projector is integral. Exact integer idempotence checks verify every
orbit projector. A verified bound on the numerical basis distance to an exact
orthonormal basis controls projected-Hamiltonian construction error. This avoids
both irrational symbolic eigenvectors and an unmanageable 73789-square Sz block.

## Outcome

The dense route with exact projectors certified δ = 9/10 at (2, 7/2) and δ = 7/10 at the odd
seed (3, 17/4), using the positive cell transfer above (proof in METHOD.md). δ = 3/5 and 1/2
have no seed within the six-cell cap; the extrapolation in METHOD.md puts their first seeds
at about 5 and 7 cells, which needs rings of 20–28 spins.
