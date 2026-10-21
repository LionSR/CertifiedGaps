# Weighted Poisson history generalization

## Candidate routes

1. Transport a general time measure through a quantile map to Lebesgue time. This gives useful
   distributional identities but makes ordered-history equality and support bookkeeping indirect.
2. Prove chamber symmetry directly for finite products of an atomless time measure. Coordinate
   collisions are null by product integration; coordinate permutations preserve product measure.
   This retains genuine time coordinates and all finite precedence/shuffle constraints.
3. Derive weighted histories only by a density relative to the existing Lebesgue construction.
   This would exclude singular continuous time laws, so it is insufficient as the main route.

## Chosen route

Use direct finite-product chamber decomposition on a measurable time domain in the real line,
then assign each finite label word the product of its nonnegative label intensities. Summing
word weights gives the power of the total intensity; summing event counts gives an exponential.
The construction should include empty domains, no labels, zero label weights, and degree zero.
Prove the Poisson count law and exponential majorants from the constructed measure rather than
assuming series convergence. Recover the existing bounded-interval Lebesgue measure by equality.
Finite-family shuffle regrouping is a subsequent stage if the generic measure API permits it.

## Cost

The proofs are uniform in event count, domain, and finite label cardinality. No numerical
spectral computation, sector decomposition, or model-specific trace identity is involved.

## Outcome

Completed five modules: generic atomless chamber decomposition, generic finite-family shuffle
integrals, weighted single-history measures, weighted finite-product absolutely convergent
shuffle expansions, and regression/axiom audits. The time-domain API is a measurable subset of
the real line. Finite-degree chamber and shuffle identities allow any sigma-finite atomless time
law. Normalization and infinite history sums use a finite time law; equality under restriction
to the allowed domain also covers a sigma-finite law with finite mass on that domain.

The full weighted product theorem includes zero-intensity words. Its proof splits the zero
product-weight case before transferring integrability to the unweighted simplex, so it does not
silently assume that a measure-zero component has an integrable unweighted observable.

Validation on October 19, 2026:

- All five modules compile with autoImplicit and relaxedAutoImplicit disabled, warnings as
  errors, and the pinned mathlib standard linter set enabled.
- The enforcing transitive axiom whitelist checked 177 declarations, permitting only propext,
  Classical.choice, and Quot.sound.
- Exact kernel replay passed for all five modules in separate processes using the unchanged
  canonical code/lean/scripts/KernelCheck.lean driver.
- Regression checks cover empty histories, empty label sets, all-zero intensities, a specific
  zero-intensity word, null domains, original Lebesgue measure equality, nonuniform disconnected
  time support, unequal weights, exact Poisson count rate, and a nonzero full shuffle expansion.
- Independent review passed: all five author modules and the independent 38-theorem regression
  module strictly compile and exactly replay. Its transitive whitelist checks 314 declarations.
  Additional tests include a nonsymmetric three-cycle, three genuine shuffles, empty-family
  integrals, the full ungrouped convergence bound, and direct finite-domain restriction laws.

The final integrated tree also passed exact kernel replay of all 111 built repository modules.
The original-to-restacked prerequisite change was verified to affect only documentation; its
Lean, dependency, and workflow blobs are identical. No shared dependency artifacts were modified.
No physical
model certificate, trace identity, or sector-specific assumption is claimed.
