# Continuous spin-one history kernel

This prerequisite constructs the continuous spin-one history kernel and its normalized complex-L²
operator. Hilbert–Schmidt bounds, compactness, and self-adjointness are proved from the actual
chronological matrix product, without additional axioms or seed hypotheses.

No certificate: this is an analytic prerequisite, not a numerical spectral-gap certificate.

## Model

A history is a finite increasing list of event times in `(0,b)`, labeled by `Fin 3`. The disjoint
union of labeled simplexes has mass `exp(3*b)`. The kernel is the real matrix trace of the merged
chronological word in `SpinOne.G`, complexified to act on `Lp ℂ 2 (historyMeasure b 3)`.

## Reflection positivity

The existing Cartesian spin-one factors are used unchanged. Exchanging the two histories relabels
their union; it neither reverses nor commutes matrix factors. This proves symmetry of the real
kernel and hence Hermitian symmetry on the complex Hilbert space.

## Seeds and certified bounds

No numerical seeds are required. The kernel has absolute value at most three, and the squared
Hilbert–Schmidt bounds are `9 * exp(6*b)` for `K b` and `9 * exp((6 - 2*a)*b)` for
`X a b = exp(-a*b) • K b`. Both operators are compact and self-adjoint for every `a : ℝ` and
`b : ℝ≥0`, including the empty-history space at `b = 0`.

## Reproduce

From the repository root with the pinned toolchain and cached dependencies:

```sh
lake build CertifiedGaps.Transfer.ContinuousKernel.Audit CertifiedGaps.Transfer.SpinOneKernelAudit
KERNEL_CHECK_MODULE=CertifiedGaps.Transfer.SpinOneKernel \
  lake env lean lean/scripts/KernelCheck.lean
KERNEL_CHECK_MODULE=CertifiedGaps.Transfer.SpinOneKernelAudit \
  lake env lean lean/scripts/KernelCheck.lean
```

The combined prerequisite modules compile in about one minute in the cloud development environment.
Kernel replay should be run once per exact module using the existing driver.

## What is rigorous

Lean checks the actual integral formulas, joint measurability, collision nullity, matrix bound,
self-adjointness, Hilbert–Schmidt estimate, and compactness. The audit enforces the standard-axiom
whitelist `propext`, `Classical.choice`, and `Quot.sound`, and tests an interleaved noncommuting word
whose reversed trace has the opposite sign.

No operator-trace, physical partition-function, spectral-moment, rotation, or sector identity is
claimed. The existing finite-dimensional gap proofs remain a separate route.

## Layout

- `IDEAS.md`: route comparison and outcome
- `../../lean/CertifiedGaps/Transfer/ContinuousKernel/`: bounded-kernel analysis and matrix bounds
- `../../lean/CertifiedGaps/Transfer/HistoryCollision.lean`: null cross-history collisions
- `../../lean/CertifiedGaps/Transfer/MergedHistory.lean`: actual chronological kernel
- `../../lean/CertifiedGaps/Transfer/SpinOneKernel.lean`: concrete operators and bounds
- `../../lean/CertifiedGaps/Transfer/SpinOneKernelAudit.lean`: regressions and axiom policy
