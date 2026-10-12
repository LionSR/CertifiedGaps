# Candidate routes

The model is the spin-1 antiferromagnet plus D(Sz)^2, with a one-site cell.
All new work stays here; the shared core is read only.

1. **Exact local RP algebra.** Integer Cartesian generators give the bond identity
   immediately. Verify the spin commutators, Casimir, onsite term and a local basis
   change, rather than assuming that the supplied generators describe the model.
2. **Symmetry-resolved dense ED.** Sz and translation reduce the largest L=12 block
   from 73789 to about 6200. Spin inversion pairs opposite Sz sectors; complex
   conjugation pairs momenta. Cost is cubic in block dimension, exponential in L;
   L=14 has blocks near 44000 and is unlikely to be feasible densely on this host.
   A floating-point scan may use low eigenvalues, but certification must use all
   levels of every block with an asserted total dimension 3^L.
3. **Verified eigensystems.** Enclose block construction, residuals and Q*Q-I.
   Bound distance to the unitary polar factor and use Weyl to enclose each sorted
   eigenvalue. Directed interval elementary functions propagate the partition sums
   and yield exact rational defect upper bounds for the core recurrence.
4. **Large-D perturbation / cluster expansion.** The product vacuum is a useful
   reference, but elementary norm perturbation is extensive and cannot certify the
   required finite-size free-energy cancellation. A convergent polymer remainder
   would need substantial new analysis, so it is not the primary route here.
5. **Smaller-seed coupling identity.** Potentially halves the largest ring but
   requires another spatial seed at 2 beta. Keep as an alternative if dense ED is
   prohibitive; do not silently substitute it for the requested thermal seed.

Chosen route: exact RP, numerical parameter screening, then complete dense spectra
and verified residual bounds for the best reachable seeds. Only even n0 can invoke
this repository's stated theorem. Numerical references and extrapolations remain
explicitly separate from certificates. A smallest successful scanned D is not a
proof that all smaller real D fail, nor a phase transition.
