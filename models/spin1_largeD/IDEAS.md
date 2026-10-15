# Candidate routes

## Follow-up plan (2026-10-14, before any new computation)

This checkout already contains the L=12 ED baseline and two certificates. Audit those
inputs before using them. The goal of this follow-up is to extract larger seeds from
the same finite inputs, with a documented failure if the bounds do not close.

* **Closed forms:** the exchange terms do not commute with the single-ion term or with
  adjacent bonds. No free-particle or frustration-free partition function is known for
  this D interval. The D=infinity product limit alone does not supply finite-D seeds.
* **Cluster/low-temperature expansion:** the product vacuum is a natural reference, but
  pair creation by exchange requires a polymer estimate uniform in volume and beta.
  An extensive remainder cannot control log Z(2n)-2 log Z(n). Deriving a convergent,
  numerically useful remainder near D=3 is not assumed. Cost after such a derivation
  could be independent of n up to exponentially small winding polymers.
* **Transfer moments (chosen beyond-ED route):** for even n, the normalized X_beta^n is
  positive. Squaring improves its purity defect s to at most
  f(s)=s^2/[2(1-s)^2]. Bound spatial purities at n=12,24,... from n=6 without constructing
  those rings. At each step use the exact identity
  P_tau(2n,b)=P_s(n,2b) P_tau(n,b)^2/P_s(n,b)^2 and P_s<=1 to bound the thermal input.
  Maintain spatial defects at both b and 2b. This costs O(log(n/6)) scalar interval
  operations after evaluating the stored L=6,12 partition functions. If the thermal
  loss prevents admissibility, record that loss explicitly, rather than calling ED
  alone the new result. A Birkhoff contraction estimate is less direct because no
  strictly positive finite kernel with controlled Trotter error is supplied here.
* **Perturbation (second chosen route):** V=sum (Sz)^2 has 0<=V<=L. Centering at L/2
  gives norm(V-L/2)<=L/2. Weyl therefore transfers each verified finite-ring spectrum
  to nearby D, with radius enlarged by L*abs(delta D)/2 and a cancelling energy shift.
  This can certify parameter intervals and a refined grid point without another ED
  run. Its defect loss is O(n beta abs(delta D)), so it will not be presented as a
  uniform-volume perturbation theorem. The subsequent bootstrap is uniform in L.
* **Smaller seeds:** the coupling identity above supplies the missing thermal purity
  on the doubled ring from smaller rings; use exact identities and track all errors.
  It may permit seed n0>6 even though ED stops at L=12.
* **Symmetry/counting:** retain Sz, translation, spin flip and paired momenta. Dense
  L=12 is a verified baseline; L=14 dimensions will be checked before any allocation.
  Because n0 must be even, L=14 alone cannot certify n0=7. Sparse numerical L=14
  diagnostics are useful only if they affect a decision. Full spectra scale
  exponentially, and unverified Krylov tails cannot certify a seed.
* **Variational sandwiches:** a trial MPS gives an upper E0(2n), but the existing
  finite spectra do not give a matching lower E0(n) with exponentially small surface
  error. A norm/Anderson bound with O(1) error will not resolve the spatial defect
  near threshold. Developing an independently verified MPS/cluster sandwich is left
  as a separate route, not silently approximated by ED extrapolation.

First replay the existing checker in this directory's own Julia project. Then test
the transfer-moment lift and centered Weyl transport, preserve exact strings in new
certificates, record source hashes and runtime/package metadata, and add focused tests.
Retain the quarter-step numerical scan as explicitly floating-point evidence; refine
the certified threshold without claiming a global minimum over D or all seed sizes.
All individual runs have a 2700-second timeout. No shared-core edits are planned.

## Original ED baseline plan

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
