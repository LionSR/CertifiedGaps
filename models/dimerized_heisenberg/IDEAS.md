# Candidate routes and first prerequisite

1. Dense ED in Sz, cell momentum, and spin-flip blocks, then rigorous residual and
   orthogonality enclosures: exponential in the largest ring (4 n0 spins for a
   seed at n0 two-spin cells). The requested 20-spin ceiling allows n0 <= 5,
   not n0 <= 10. Known n0=4 seeds require 16 spins.
2. Expansion around disconnected dimers: potentially scalable in n0, but a
   rigorous finite-temperature remainder is needed away from delta=1.
3. Generalized cell reflection, using swap composed with time reversal: may
   support a broader transfer theorem, but is not automatically the real-C
   identical-factor decomposition required in the task.

Chosen first step: verify the exact local RP assumption before ED. A bond of
form sum c C tensor C must commute with cell exchange. The physical oriented
inter-cell Heisenberg bond does not. Common cell basis changes preserve this
obstruction. A sublattice rotation is a common cell unitary and cannot fix it.
If confirmed, stop certification at hypothesis (a); numerical seeds cannot
repair an algebraically impossible identity. Do not label this as a failure to
find an admissible numerical seed or as evidence against a physical gap.

## Follow-up scope correction and numerical route

Independent review confirms the general repository criterion IS available:
Theta is internal cell reflection composed with conjugation. The outer-cell
swap obstruction excludes only the stronger literal identical-factor demand.
Larger even blocks preserve this reflected form; alternating orientations give
two different identical-factor bond operators, not one homogeneous operator.

Implemented full dense Sz/momentum/spin-flip ED through 16 spins for all six
parameters, with complete dimension accounting and independent small-ring
checks. Scanned n0=2,3,4 and beta=0.25:0.25:32; odd n0 remains diagnostic only.
This supplies useful non-certified data despite the original specification
mismatch. Residual enclosures and rational purity certification remain undone.

## Outcome

Verified enclosure route completed.

Following explicit authorization to continue under an audited applicable
repository theorem, implement interval algebraic block construction and dense
LAPACK eigenvectors. Compute residual and orthogonality centers using four
real GEMMs; bound every rounding/underflow contribution by exact rational
gamma bounds. A polar-factor/Weyl theorem then bounds every true eigenvalue
with multiplicity. Outward full-spectrum log partition sums give exact
rational purity seeds and positive general-reflected-RP bounds for delta
0.5, 0.3, 0.25. The literal original real identical-factor requirement remains
false and explicitly excluded from the claim. Independent review validated
the asymmetric recurrence/tail directly from the repository note identities.
