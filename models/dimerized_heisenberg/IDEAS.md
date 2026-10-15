# Current task: route assessment before computation (2026-10-14)

The existing results below are historical evidence, not trusted inputs for this task.
Re-derive the conjugated-factor RP identity, rerun verified spectra, and make the main checker
recompute spectral enclosures rather than trust saved residual maxima. The required form is
`-h = sum C tensor conj(C)`: complex channels are allowed. The old stronger real-channel
obstruction is irrelevant to the current task and must not be presented as a failed hypothesis.

The computation ceiling is **10 cells = 20 spins in the largest ring**, not a seed of 10
cells requiring 40 spins. Direct seeds have even n0, hence n0=2 or 4. An n0=5 diagnostic can
use 10 cells but cannot establish the stated even-seed theorem. All runs get a 44-minute
external timeout. No existing floating-point survey value will be used as a proof input.

Candidate routes, including ones that can exceed the ED ceiling:

1. **Closed forms.** At delta=1, independent strong dimers give
   `Z_n=(exp(3 beta/2)+3 exp(-beta/2))^n` at all sizes in O(1) operations. At the requested
   deltas the alternating isotropic Heisenberg chain is neither commuting nor free fermionic;
   Jordan-Wigner retains the density interaction. No closed form is assumed.
2. **Rigorous expansions.** A linked polymer expansion around dimers would bound the
   finite-ring remainder uniformly in n. Its expansion parameter involves the weak/strong
   ratio `(1-delta)/(1+delta)` (1/3 at delta=1/2, 7/13 at delta=3/10, and larger below),
   imaginary-time integrals and connectivity. High-temperature convergence conflicts with
   the low temperature needed for thermal purity. Deriving a valid low-temperature cluster
   remainder is substantial work; merely truncating a series is not certification.
3. **Transfer operator.** A rigorous subleading/leading eigenvalue ratio would make the
   spatial defect exponentially small at arbitrary n, with cost controlled by the transfer
   approximation and beta rather than Hilbert dimension. RP supplies self-adjointness but
   no contraction constant. A Trotter/MPO estimate needs both discretization and truncation
   bounds; a positive-kernel/Birkhoff estimate is not supplied by the spin interaction alone.
4. **Perturbation from dimers.** The exact inequality `|log Z-log Z_dimer| <= beta ||V||`
   with `||V|| <= 3 n (1-delta)/4` scales as beta*n*(1-delta). Combining independent ring
   errors loses cancellation of the bulk free energy, so increasing n makes this elementary
   bound worse. Record its certified failure at requested seeds; it still gives a scalable
   sufficient condition close to delta=1. A linked remainder could fix the volume loss.
5. **Smaller-seed identities.** Evaluate the precise hypotheses before using the note's
   thermal replacement. Bounds on `P_tau(n,beta)` and `P_s(n,2beta)` still involve a doubled
   ring unless the latter is independently bounded by a transfer/expansion argument. They
   cannot silently license an odd seed or a 10-cell seed within a 20-spin ED ceiling.
6. **Symmetry/counting.** Use Sz, translation by one cell and spin-flip orbits exactly.
   Their dimensions and multiplicities must exhaust `2^(2n)`. Complete dense spectra cost
   roughly the sum of block-dimension cubes, with memory quadratic in the largest block.
   SU(2) character inversion or exact trace recursions could go further, but require a
   remainder bound for uncomputed energies/traces. Partial Lanczos alone cannot bound Z.
7. **Variational sandwiches.** Trial states upper-bound E0(2n); cluster/Anderson bounds
   lower-bound E0(n). For spatial purity a useful estimate must also control excited-state
   partition tails; independent extensive energy errors again lose the cancellation needed
   at beta of order 1/gap. No numerical MPS energy is accepted as a lower bound.

Chosen route: first repair the proof/checker boundary and regenerate complete verified
8/16-spin spectra (ED is the baseline); accompany it with an exact, scalable dimer-limit
perturbation bound and documented assessment of why it does not reach the requested deltas.
Then scan all specified deltas with complete symmetry counting and extend diagnostics toward
20 spins within each time cap. Negative findings must distinguish a tested beta grid, an
all-beta exclusion if one is proved, and resource-limited cases. Ratios to a finite-ring gap
must be labeled accordingly; no finite-size value is an exact thermodynamic Delta. Compare
seed sizes to delta^(-2/3) only as a scaling diagnostic, not a rigorous extrapolation.

## Historical memo (retained to explain the existing implementation)

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
