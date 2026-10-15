# Task 5: small seeds, exact arithmetic, and the remaining trust boundary

The cheapest convincing route found is the longitudinal Ising point
`g=1/4, h=1/20, n₀=6, β₀=97/4`. Two approximate eigenvectors per ring, together with a
classical diagonal comparison, suffice numerically. Full Rayleigh–Ritz matrices, a
recursive half-chain anchor, momentum sectors, and SU(2) are unnecessary for this point.

**This deliverable does not close an unconditional RPRing seed.** It supplies actual
kernel-checked arithmetic, the exact original seed conclusions under explicit partition
enclosure hypotheses, and a quantitative feasible route between them. No original model
theorem is newly unconditional. An axiom audit of a conditional theorem does not discharge
its hypotheses.

All new work is under `models/lean_seed_small/`. The original `lean/` and
`models/spin1/results` were not modified. The project uses Lean `v4.35.0-rc3` and Mathlib
commit `c55e6e786f49471c72fbddbec5415808896aec1e`, matching the original manifest.

## What was actually checked

`certify.jl` uses Julia 1.12.1/LAPACK only to search for vectors. It rounds them to integers,
with scale `S=2^52`, and emits exact numbers as strings and Lean source. The certificate is
about 142 kB of JSON; its generated Lean constants are about 141 kB. These are vectors,
not a claimed kilobyte-size Krylov recipe. The padded packed representation is reconstructed
from compact 56-bit slots by twelve masked shifts; padding is not stored.

For the integer Hamiltonian `A=20H`, the checked quantities are

```
|vᵢᵀvᵢ − S²| ≤ 10⁻¹² S²
|v₁ᵀv₂|       ≤ 10⁻¹² S²
Σᵢ ‖S A vᵢ − λᵢ vᵢ‖² ≤ 10⁻¹⁸ S⁴
λ₁ + S < λ₂,       λ₂ + S < (20ρ) S.
```

The physical candidate energies are `μᵢ=λᵢ/(20S)`. The last two tests provide margins of
more than `1/20` for distinctness and separation from the anchor.

| Ring | Dimension | Candidate E₀ | Candidate E₁ | Proposed anchor ρ | Exact-residual norm / (20 S²), displayed approximately |
|---|---:|---:|---:|---:|---:|
| 6 sites | 64 | −6.391774806847205 | −5.796599647006008 | −3.7 | 5.172×10⁻¹⁵ |
| 12 sites | 4096 | −12.78354960854820 | −11.59319932209873 | −11.5 | 7.475×10⁻¹⁴ |

The decimal display and LAPACK eigenvalues are **non-certifying**. The integer inequalities
in the Lean theorems are certifying statements about the explicitly defined arithmetic.

- `Generated.check6` and `check12` evaluate the packed arithmetic using `decide +kernel`.
- `Generated.direct6` evaluates the scalar Hamiltonian recipe directly.
- `Generated.anchor6` and `anchor12` count classical diagonal entries by direct enumeration.
- `GeneratedChunks` checks the 4096-state scalar calculation in 64 separate chunks, avoiding
  the excessive peak memory of one scalar reduction. It also checks the enumeration coverage.
  `ChunkCoverage` proves that the generated chunk list is exactly the fixed range enumeration;
  a generator cannot silently duplicate or omit chunks. This gives a direct arithmetic route
  independent of proving the general packed backend correct.
- `Arithmetic.spatial_arithmetic` and `thermal_arithmetic` check the exact rational endpoint
  inequalities used by the original seed targets.
- `ResidualLemma` proves the finite weighted-sum step: a residual bound forces at least one
  spectral value into the specified interval. Matrix diagonalization/norm identities are not
  supplied by that lemma.
- `SeedBridge` concludes exactly
  `1−9064e−9 ≤ (IsingLongitudinal.model (1/4) (1/20)).spatialPurity 6 (97/4)` and
  `1−3141e−8 ≤ (IsingLongitudinal.model (1/4) (1/20)).thermalPurity 12 (97/4)`,
  **assuming** its four-field `LogEnclosures` structure.

The imported original module contains seed sorries, but these new results do not invoke
them. The axiom outputs, build status, and measured timings are recorded in `VALIDATION.md`.
The full Lake project passed its final build; all 17 audited declarations use only the
standard axioms. `first_point_seed_instance` connects the displayed numeric goals to the
literal first entry of the original `points` list. [AXIOMS.txt](AXIOMS.txt) preserves the
actual axiom output. Its enclosure hypothesis remains open.

## A strictly cheaper route for the selected point

Let `D` be the classical diagonal part of `H`. Since `−gX ≥ −gI` on one spin,
`H ≥ D − ng I`. There are exactly two exceptional classical configurations: all up and all
down. Every other periodic configuration has at least two domain walls. Thus its classical
energy is at least `−n(1+h)+4+2h`. At `g=1/4, h=1/20`, this gives

```
ρₙ = −n(1+h+g)+4+2h = (82−26n)/20,
λ₃(H₆) ≥ −3.7,          λ₃(H₁₂) ≥ −11.5.
```

The direct classical counts are proved in Lean. The operator inequality and eigenvalue-count
consequence are still mathematical proof obligations. The anchor is not accepted merely
because the diagonal count passed.

An individual residual gives an eigenvalue within
`10⁻⁹/(20√(1−10⁻¹²)) < 5.01×10⁻¹¹` of each candidate, once the norm identities are proved.
The two intervals are disjoint and below ρ. A proved count of at most two eigenvalues below
ρ would then make them exactly the first two. **Cross-Gram orthogonality is unnecessary** for
these disjoint intervals; it is checked here for comparison with the proposed plan.

This removes the full polar/Weyl Gram-matrix argument. For larger spectra the analogous
optimization is to certify small clusters of overlapping eigenvalue intervals and their
multiplicities, rather than form a Gram matrix across all k levels. That optimization still
needs a counting theorem and a reliable anchor; it is not a license to discard uncounted states.

With a deliberately conservative energy radius `10⁻¹⁰`, the two-level trace sandwich is

```
Σᵢ₌₁² exp(−β(μᵢ+r)) ≤ Zₙ(β)
Zₙ(β) ≤ Σᵢ₌₁² exp(−β(μᵢ−r)) + (2ⁿ−2) exp(−βρₙ).
```

The following is a **non-certifying Float64 feasibility experiment**, reproducible with
`feasibility.jl`; it assumes the spectral theorem described above.

| Quantity | Sandwich result | Required target |
|---|---:|---:|
| Spatial defect upper bound | 1.210666×10⁻⁶ | 9.064×10⁻⁶ |
| Thermal defect upper bound | 9.949077×10⁻⁹ | 3.141×10⁻⁵ |
| 6-site normalized tail at β | 2.7771×10⁻²⁷ | included above |
| 12-site normalized tail at β | 1.2425×10⁻¹⁰ | included above |

The predicted log-partition intervals also fit strictly inside the four existing external
endpoints used by `SeedBridge`. Therefore the remaining task is a proof task, not a request
for sharper numerical vectors. The `h=1/10` point has only one level below the inventory
cutoff, but the simple rank-two classical anchor does not by itself isolate that one level;
the `h=1/20` point is the cheapest route found with a complete numerical anchor argument.

A concrete scalar-proof route is `Real.exp_bound` on each exponent divided by 1024
(the largest exponent here is about 620), followed by ten outward-rounded rational
squarings. A 16-term Taylor proposal is amply plausible at the required margins, but its
rounding and remainder checks must themselves be proved. This evaluates rational arithmetic
in Lean, not a numerical implementation of real `exp`.

## Inventory: dimensions, precision, and spectrum sizes

[TABLES.md](TABLES.md) lists **every one of the 43 parameter points**, including both epsilons,
β, dimension, tail budget, planning eigenvalue radius, cutoff, level count, and arithmetic cost.
[inventory.tsv](inventory.tsv), [inventory.json](inventory.json), and [costs.json](costs.json)
contain the reproducible data and source hashes.

The seed defects `εₛ, εₜ` are not eigenvalue tolerances. To avoid silently using all of a seed's
slack, the inventory sets
`σ=min(εₛ−εₛ,source, εₜ−εₜ,source)`, `η=min(10⁻⁶,σ/16)`,
`r=η/(4β)`, and counts levels below `E₀+(log D+log(1/η))/β`.
Here `source` means the stored external enclosure, **not a trusted Lean theorem**. This is a
conservative planning convention, not a necessary precision theorem. The selected point
demonstrates why recomputing tighter bounds can permit a much larger error budget than this
source-slack convention. Counts use the nominal stored spectra with multiplicities and are
non-certifying; no external radius or BLAS error model is trusted by Lean.

| Family | Largest full dimension | Largest Sz block, when used | k at the stated cutoff | Full three-Gram arithmetic scenario |
|---|---:|---:|---:|---:|
| Longitudinal Ising, 30 points | 4,096 | 4,096 (no Sz reduction) | 1–145 | 0.014–150 seconds per point |
| Spin-one dimerized δ=.9 | 6,561 | 1,107 | 195 | 5–80 seconds |
| Staggered Heisenberg, 5 points | 65,536 | 12,870 | 424–1,249 | 0.14–19.8 hours |
| Heisenberg ladder J=4 | 65,536 | 12,870 | 1,238 | 1.16–19.4 hours |
| Spin-half dimerized, 3 points | 65,536 | 12,870 | 862–1,558 | 0.56–30.8 hours |
| Ising ladder g=3 | 65,536 | 65,536 (Sz is not conserved) | 322 | 15.7 hours |
| Large D=4 | 531,441 | 73,789 | 1,480 | 15.4–379 hours |
| Large D=3.206 | 531,441 | 73,789 | 2,465 | 42.5–1,051 hours |

These are **arithmetic-only scenario estimates**, not measured certificate runtimes or proven
lower bounds. For conserved charge they range from equally distributed k across charge sectors
to all k concentrated in the largest block, pricing each product at the largest block width.
They exclude vector generation, orthogonalization, packing, anchor construction, spectral
proofs, and scalar exponential bounds. Smaller actual blocks or clustered certificates can
reduce them substantially. No momentum, cyclotomic, or SU(2) saving is assumed.

The low β of the large-D seeds is why thousands of levels matter despite their apparently
small seed lengths. These results reject a blanket 15–60 minute estimate for these models
with full k×k Gram checks. They do **not** by themselves refute the separate spin-one Heisenberg
proposal at β=49/4, whose low-energy count is different.

## Measurements and changes to the implementation plan

The primitive benchmark uses two distinct dense Nat operands `(2ᵇ−1)(2ᵇ−3)`, checked against
`2²ᵇ−2ᵇ⁺²+3` by `decide +kernel`. The measured kernel type-checking times were:

| Operand bits | Kernel time |
|---:|---:|
| 262,144 | 0.00111 s |
| 1,179,648 | 0.00471 s |
| 3,706,560 | 0.0304 s |
| 21,251,232 | 0.415 s |

They are primitive microbenchmarks, with no claim that arbitrary certificate checks attain
this throughput. `costs.jl` uses a local size exponent of 1.15 and the nearest measured size.

The direct single-expression 12-site checker was stopped after reaching 7.4 GB and 38 CPU
seconds without completing. The first packed version checked the large vector arithmetic
in about 0.4 seconds but spent about 140 seconds elaborating two padded decimal constants.
The compact expansion recipe cut that elaboration to about 6 seconds and the packed check
to about 0.3 seconds in a later run. Direct classical counting took about 52 seconds.
Separate scalar chunks provide another route for the small target without requiring a generic
packed correctness library. See `VALIDATION.md` for their final build measurement.

Changes to the proposed plan:

1. Use the two-state residual/counting route for the first Ising point. Do not implement a
   recursive anchor or full Ritz/Gram machinery for it.
2. Make scalar exponential enclosures an explicit proof component. Julia can propose rational
   endpoints, but an endpoint table alone does not certify `exp`, `log Z`, or a spectral sum.
3. Include a proved bridge from the computational-basis integer matrix to `RPRing.hamiltonian`,
   and from executable sums to norms. Exact residual arithmetic alone does not supply either.
4. Store compact coefficients and expand padding in Lean. Literal elaboration can cost far
   more than the GMP multiplication being benchmarked.
5. Do not treat packed primitive timings as full-checker timings. Count vector preparation,
   all Gram products, classical counts, scalar bounds, and peak kernel memory separately.
6. For large k, investigate disjoint clusters before implementing all k² Gram products.

## Remaining proof work and revised estimates

The exact missing theorem for the selected original model is `LogEnclosures`. Its fields are
visible in `SeedBridge.lean`; none has been assumed as a new global axiom. A route to proving it
has the following remaining pieces:

| Component | Additional Lean LoC estimate | Concrete obligation |
|---|---:|---|
| Matrix identification and scalar-sum interpretation | 200–400 | Binary states ↔ `Fin n → Fin 2`; prove `A=20H`; identify checked norms/Gram entries |
| Residual-to-spectrum bridge | 150–350 | Apply the proved weighted residual lemma through the Hermitian spectral theorem |
| Classical comparison and counting | 250–500 | Prove `H≥D−ngI`, transfer the checked classical count via min–max/inertia, and exhaust two eigenvalues |
| Rational exponential/trace bounds | 200–450 | Prove outward scalar exp bounds and the finite trace sandwich; derive the four log bounds |
| Integration and proof checks | 100–250 | Convert exact margins, connect to the original chain without seed sorries |
| Optional reusable packed backend proof | 300–700 | Mask/flip/reverse, carry bounds, middle-product dot extraction, widening; not needed if the direct chunk route is used |

Thus the next complete selected-point proof is estimated at **900–1,950 additional Lean lines**
using the direct route, or **1,200–2,650** with a reusable packed backend. This is a task breakdown,
not a measured completion guarantee. A reasonable planning range is 1–2 weeks for an experienced
Lean contributor, with the counting/min–max bridge the main uncertainty. Once those lemmas exist,
a 1–5 minute kernel run for this particular point is plausible from the measured arithmetic;
the complete proof has not been timed.

For a reusable backend covering all seven model families, 4.5–7.5k lines is optimistic: explicit
basis/symmetry identification, anchors, cluster multiplicities, and scalar bounds need a separate
budget. A planning allowance of **8–15k additional lines and 4–8 weeks** is more defensible for
that broader scope, conditional on finding usable anchors. The larger Gram scenarios rule out
a uniform minute-scale runtime without further structural reductions.

## Which models become unconditional, and in what order

**None became unconditional in this task.** All listed targets are plausible finite certification
problems; this work establishes feasibility evidence rather than an unconditional existence claim
for their proposed efficient certificates.

1. The selected longitudinal Ising point. Its explicit two-state anchor and large seed margin
   make it the best immediate formalization target. A separate theorem for this point can be
   unconditional before the universal 30-point `seeds` theorem is finished.
2. The remaining longitudinal Ising points and spin-one dimerization δ=.9. The Ising family
   reuses matrix identification; δ=.9 has a small 1,107-dimensional largest Sz block. General
   points still need anchors. Proving an Sz-basis transformation is extra work for spin one.
3. Strong staggered fields (h=1, 3/4), then the Heisenberg ladder and the easier spin-half
   dimerization targets. Reuse the spin-half matrix/symmetry infrastructure and prioritize by
   actual cluster counts and anchor strength, not full dimension alone.
4. Weaker staggered fields, harder dimerization targets, and the Ising ladder. The latter has
   no Sz decomposition and its full Gram scenario is particularly expensive.
5. Large D=4, then D=3.206. They have the largest level counts; D=3.206 also has only about
   `5.85×10⁻¹⁰` slack relative to its stored enclosure. Redesigning the seed or using clusters
   is more promising than blindly applying the dense Gram plan.

Open risks are the unproved matrix/sum bridge, min–max/counting formalization, analytic scalar
enclosures, symmetry multiplicities for other models, anchor quality, and runtime growth from
recipe regeneration. The generator's correctness and BLAS rounding model are not soundness
assumptions: wrong proposals must merely fail the kernel checks.

Final validation: the Lake build, 17-declaration axiom audit, 8 Julia tests, pinned versions,
generator hash, and protected-source checks passed. The delivered outcome is a reproducible
partial proof package, not an unconditional gap theorem.
