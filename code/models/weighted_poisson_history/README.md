# Weighted Poisson histories

This analytic prerequisite constructs weighted chronological event histories for arbitrary
finite atomless time measures on the real line. It proves the factorial chamber law, exponential
mass, normalized Poisson event-count distribution, and absolutely convergent finite-family
chronological shuffle expansion from the actual measures.

**No certificate: analytic prerequisite.** No Hamiltonian, spectral gap, sector decomposition,
or partition-function trace identity is asserted.

## Model

An arbitrary finite label type carries nonnegative intensities `w`. An arbitrary finite atomless
Borel measure `ν` on the real line supplies event times in a measurable domain `s`. The domain
may be the entire real line, disconnected, or empty; the time law need not have a density.

A history consists of a finite label word and a strictly increasing time tuple in `s`. Its
component measure is the product label intensity times the restriction of the actual product
time measure to that ordered simplex. Empty histories are included. Labels of zero intensity
and empty label types are allowed.

The underlying chamber and finite-degree shuffle theorems are stronger: they require only
sigma-finite atomless `ν`, with no finiteness assumption on `ν(s)`.

## Reflection positivity

Not applicable. The observable can be Banach-valued, and chronological reorganization never
commutes or exchanges its values. Reflection positivity and matrix trace identities remain
separate obligations for any later physical application.

## Result

Let `Λ` be the sum of the label intensities and `T = ν(s)`.

- Each ordered degree-`n` simplex has mass `T^n / n!`.
- The actual weighted history measure has total mass `exp(Λ T)`.
- Normalizing by `exp(-Λ T)` gives a probability measure whose event count is exactly Poisson
  with rate `Λ T`.
- The unnormalized event-count moment is `∫ R^N = exp(Λ T R)` for every nonnegative `R`.
- An almost-everywhere strongly measurable observable bounded in norm by `C R^N` is
  integrable and has L¹ norm at most `C exp(Λ T R)`.
- For a finite family sharing the time law and domain, with independently chosen finite label
  alphabets and intensities, the actual product history measure decomposes into count-vector,
  label-word, and admissible chronological-shuffle terms. Both the individual norm integrals
  and the norms of the individual Banach-valued terms are summable. A product exponential
  bound supplies the convergence hypothesis and its explicit quantitative majorant.
- A zero-weight word contributes zero even when its unweighted observable is not integrable.
  There is no positivity assumption silently excluding zero-intensity labels.
- Unit label intensities and Lebesgue time restricted to `(0,b)` recover both existing history
  measures, unnormalized and normalized, by actual measure equality on definitionally equal
  history spaces. Restricting any time law to its allowed domain also preserves the history
  measure by equality.

For a sigma-finite time law that is only finite on the chosen domain, apply the normalization
and product APIs to `ν.restrict s`. Its atomlessness is inherited, and its total mass is `ν(s)`.
The theorem `historyMeasure_restrict` identifies this construction with the original history
measure on `s`; applying that equality to each factor gives the identical product history
measure. Thus the global finite-measure instance in the normalization/product APIs does not
require a globally finite original law when the allowed domain has finite mass.

The regression audit uses a nonuniform time law supported on two disjoint intervals, with
masses one and two, and label intensities zero, two, and three. It checks ordered-simplex mass
`9/2` at degree three, weighted degree-two mass `225/2`, normalized count rate `15`, degenerate
cases, original-measure specialization, and a nonzero full weighted shuffle expansion.

The independent audit also checks infinite-mass cases, direct finite-domain normalization,
inverse permutations with a nonsymmetric integrand, a fixed weighted `(1,0,2)` shuffle sum
of `324`, exactly three interleavings, empty families, arbitrary-observable zero-weight
components, and an explicit ungrouped norm bound.

## Reproduce

Use the repository's pinned Lean 4.36.0-rc1 and mathlib revision
`807d3836b165983210ee2ecf3cd28f239436949c`. From the repository root:

```sh
lake build CertifiedGaps.Transfer.WeightedPoissonHistoryAudit \
  CertifiedGaps.Transfer.WeightedHistoryAudit
for module in AtomlessTimeOrdering AtomlessHistoryShuffle WeightedPoissonHistory \
  WeightedProductPoissonHistory WeightedPoissonHistoryAudit WeightedHistoryAudit; do
  KERNEL_CHECK_MODULE=CertifiedGaps.Transfer.$module \
    lake env lean code/lean/scripts/KernelCheck.lean
done
```

With cached dependencies and one Lean thread, direct strict compilation of the six new
modules took 34 seconds in the development environment (excluding dependency download/build).
The development validation additionally passes `-DautoImplicit=false`,
`-DrelaxedAutoImplicit=false`, `-DwarningAsError=true`, and
`-Dweak.linter.mathlibStandardSet=true` to every module, with `LEAN_NUM_THREADS=1`.
It uses private output artifacts and the pinned read-only dependency cache; no shared
artifacts or dependency sources are modified.

## What is rigorous

The theorem statements and proof terms are checked by Lean. No `sorry`, new axioms, numerical
approximation, `native_decide`, or assumed series convergence is used. The audit enforces a
transitive whitelist allowing only `propext`, `Classical.choice`, and `Quot.sound` across all new
namespaces. Exact kernel replay runs the unchanged bundled Lean replay routine for each new
module in a separate process. The author audit enforces the whitelist over 177 declarations;
the independent regression audit checks 314 declarations, including its 38 named regressions.

These results concern real ordered event times and finite label alphabets. They do not assert
a history-space equivalence for arbitrary other ordered measurable spaces, infinite-intensity
normalization, infinite label alphabets, trace identities, or a spectral-gap certificate.

## Layout

- `IDEAS.md`: route selection, scaling, and final validation outcome
- `code/lean/CertifiedGaps/Transfer/AtomlessTimeOrdering.lean`: product collisions, chamber mass, and
  arbitrary-integrand chronological cube decomposition
- `code/lean/CertifiedGaps/Transfer/AtomlessHistoryShuffle.lean`: precedence constraints and actual
  finite-family product-time shuffle integrals
- `code/lean/CertifiedGaps/Transfer/WeightedPoissonHistory.lean`: weighted measure, normalization,
  Poisson law, moments, domination, and actual Lebesgue specialization
- `code/lean/CertifiedGaps/Transfer/WeightedProductPoissonHistory.lean`: finite-product measure
  decomposition, weighted shuffle terms, and absolute convergence
- `code/lean/CertifiedGaps/Transfer/WeightedPoissonHistoryAudit.lean`: regression checks and enforcing
  transitive axiom whitelist
- `code/lean/CertifiedGaps/Transfer/WeightedHistoryAudit.lean`: independent numerical, degeneracy,
  shuffle, finite-domain restriction, and transitive-axiom regression proofs
