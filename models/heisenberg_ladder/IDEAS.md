# Heisenberg ladder: pre-computation route memo

Written before any model computation. Lengths below count rungs (two spin-1/2 sites).
The supplied survey values are hypotheses to reproduce, never certified inputs.

## Candidate routes and their costs

1. **Closed forms / independent rungs.** At zero leg coupling the rung singlet and
   three triplets give Z_n = (exp(3 beta J/4) + 3 exp(-beta J/4))^n.
   The interacting isotropic ladder is neither commuting nor free fermionic;
   Jordan-Wigner introduces interactions. The independent-rung formula is an
   exact arbitrary-length reference, but cannot replace the interacting trace.
2. **Rigorous expansions.** A strong-rung linked-cluster / polymer expansion could
   isolate wrapping clusters and bound log Z_2n - 2 log Z_n exponentially in n.
   Its cost would depend on expansion order and beta, not exponentially on n.
   At the requested moderate J a convergence majorant and wrapping remainder
   must first be proved; none of the quoted survey numbers provides one. A
   high-temperature expansion is poorly matched to the low-temperature seed.
   Do not claim an unproved convergence radius.
3. **Transfer operator.** In the checkerboard-rotated basis all cross-cell
   coefficients have exact real RP factors. A positive Trotter transfer kernel
   has dimension growing exponentially with imaginary-time slices. Certified
   contraction plus a Trotter remainder could reach arbitrary n, but the
   low-temperature contraction estimate may deteriorate with slice count.
   RP alone gives self-adjointness, not a quantitative contraction coefficient.
4. **Perturbation from independent rungs (selected analytical test).** Each of
   the 2n leg bonds has norm 3/4, so ||V_n|| <= 3n/2. Weyl gives the exact
   sandwich |log Z_n - log Z_n^(0)| <= 3 beta n/2. Apply this bound to both
   purities, retaining the extensive errors, and determine analytically whether
   it can close a seed for any n and beta. Cost O(1), arbitrary n. This route
   may fail because it discards the cancellation between bulk free energies.
   Record that failure explicitly; do not quietly treat an ED baseline as an
   arbitrary-size seed method. A second use of perturbation is to continue a
   verified finite seed over a coupling interval with a proven norm remainder.
5. **Smaller seeds.** Check the precise thermal/spatial inequality in the notes
   before using survey Theorem 2.24(v). Replacing one purity is useful only if
   every remaining trace and omitted sector is bounded. No odd n0 is allowed
   by the stated certificate standard; ten-rung rings are scan data, not an
   automatic five-rung seed.
6. **Symmetry and exact counting (selected finite baseline).** Work spectrally
   in the unrotated SU(2)-invariant basis with total Sz and rung momentum, and
   prove its unitary equivalence to the exact rotated RP model. Spin inversion
   and conjugate momenta provide exact multiplicities. Full verified spectra
   for four/eight rungs use the core residual/Weyl enclosure, not raw LAPACK.
   For ten rungs use sparse sector-resolved low-energy calculations with
   explicit missing-state counts; label unverified results floating point.
   Full dense cost is sum O(block_dimension^3), memory O(max block_dimension^2).
   Translation reduces the largest block by approximately n, but remains ED.
7. **Variational sandwiches.** Rung-singlet products give a rigorous trial
   energy; finite-cluster lower bounds can improve the ground-energy sandwich.
   MPS and Lanczos alone do not control the thermal tail or prove the spatial
   defect. Their cost can be polynomial in n at fixed bond dimension, but a
   small bulk energy-density error still grows extensively in the purity ratio.

## Chosen work and stopping rules

First document the arbitrary-size perturbation test. Then reproduce J=4 with
exact RP factors, complete sector multiplicities, rigorously enclosed log Z,
and the core exact rational recurrence. Scan J in {7/2,3,5/2,2,3/2} on even rings
through ten rungs. Store reproducible evidence, distinguish complete spectra
from partial low-energy data, and certify every claimed seed independently.
Use a timeout of at most 2700 seconds for each run. No floating-point scan is a
proof of nonexistence of a seed or of a gapless model.

Unless clarified otherwise, “threshold” means the threshold reached by this
finite-size certificate method, not the physical gap-closing coupling. Report
an achieved coupling and a bracket/search resolution only to the extent the
calculation supports them. A failed finite scan cannot establish the physical
threshold. All implementation and outputs stay in this model directory.

## Checkout discovery and route update (before fresh computation)

The checkout already includes a complete earlier ED implementation, stored
sector-resolved scans through ten rungs, and a J=4 certificate. The prior memo
is preserved as `IDEAS.previous.md`, including its memory observations. Reuse
that evidence explicitly as historical floating-point data; do not repeat six
long dense scans merely to overwrite equivalent data. Rebuild symmetry counts
and replay every stored scan, while freshly enclosing all spectra used by new
certificates. Strengthen the result with certified lower-coupling seeds, a
rigorous finite-beta-grid failure test, and local coupling continuation.

The centered independent-rung sandwich is stronger than the uncentered memo
estimate but still fails for every n,beta when 0<J<=4; `PERTURBATION.md` gives
a proof without numerical input. This settles the selected arbitrary-length
route negatively rather than leaving it as an untested suggestion.

Before running the new suite, add a continuous-temperature negative test for
J=2 and J=3/2: shift every verified spectrum below its ground energy, use
monotonicity of the shifted partition sums and thermal purity to bound defects
on beta intervals, and split intervals until a lower-defect recurrence proves
rejection. A positive enclosed ground-energy mismatch supplies the beta-to-
infinity tail. This makes no assumption that spatial purity is monotone.
