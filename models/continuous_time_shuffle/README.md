# Chronological history shuffles

This analytic prerequisite proves that products of ordered bond-history time simplexes
reorganize as sums over the global chronological interleavings preserving every bond's
internal order. It applies to any finite event counts and any Banach-valued integrand.

No certificate: this is a model-independent measure-theoretic result, with no numerical seed
or spectral-gap claim.

## Model

A finite family of histories has event counts `counts : ι → ℕ`, where `ι` is any finite type.
Each history has strictly increasing times in `(0, b)`. The total event degree is
`∑ i, counts i`. Zero counts and an empty family are included.

## Reflection positivity

No reflection-positivity or commutativity hypothesis is required. The value space is an
arbitrary real Banach space, so ordered matrix and operator expressions are admissible
integrands when they satisfy the stated integrability hypothesis.

## Seeds and certified bounds

There are no numerical seeds. The result is an exact identity of Lebesgue integrals:
the product of the separate simplexes equals the finite sum over admissible interleavings
of one global simplex. The inverse shuffle recovers each original event slot.

## Reproduce

From the repository root:

```sh
lake build CertifiedGaps.Transfer.HistoryShuffleAudit \
  CertifiedGaps.Transfer.ProductPoissonHistoryAudit
KERNEL_CHECK_MODULE=CertifiedGaps.Transfer.HistoryShuffle \
  lake env lean lean/scripts/KernelCheck.lean
KERNEL_CHECK_MODULE=CertifiedGaps.Transfer.HistoryShuffleAudit \
  lake env lean lean/scripts/KernelCheck.lean
```

The finite-degree source builds in approximately three seconds with the pinned Lean 4.36.0-rc1
toolchain and cached dependencies. The product-history module and its audit add a few seconds.
Use the same replay command for `CertifiedGaps.Transfer.ProductPoissonHistory` and its audit.

## What is rigorous

Lean verifies the null-collision chamber partition, the exact within-history constraints,
the finite-family product-Lebesgue change of coordinates, and the Bochner integral identity.
The audit allows only the standard axioms `propext`, `Classical.choice`, and `Quot.sound`.

`HistoryShuffle` contains the fixed finite-degree identities. The companion
`Transfer/ProductPoissonHistory.lean` now proves the actual finite-product history measure
decomposition and the fully ungrouped, absolutely convergent count/label/shuffle series, with
an explicit exponential majorant. Constructing the kernel operator and identifying a physical
partition function are separate components.

## Layout

- `IDEAS.md`: route comparison and outcome
- `../../lean/CertifiedGaps/Transfer/HistoryShuffle.lean`: proofs
- `../../lean/CertifiedGaps/Transfer/HistoryShuffleAudit.lean`: finite-degree regressions and audit
- `../../lean/CertifiedGaps/Transfer/ProductPoissonHistory.lean`: full product-history series
- `../../lean/CertifiedGaps/Transfer/ProductPoissonHistoryAudit.lean`: product-series regressions
