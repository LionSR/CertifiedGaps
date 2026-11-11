# Continuous spin-one history kernel

## Candidate routes

1. Construct the chronological kernel on the actual finite-history measure, prove its pointwise
   bound and Hermitian symmetry, and derive the integral-operator properties from general
   finite-measure bounded-kernel results. This isolates the analytic and combinatorial obligations.
2. Approximate the kernel by finite time meshes, then identify an operator-norm limit. This adds
   discretization and convergence obligations that are unnecessary for the bounded history kernel.
3. Assume the desired kernel is already a compact or Hilbert–Schmidt operator. This does not
   construct the paper's operator and would leave the main prerequisite unproved.

Choose route 1. Define the kernel as a finite sum of chronological chambers on each pair of history
components. At most one chamber contributes, so the trace bound has no factorial loss. Swapping
histories relabels events and leaves their ordered matrix word intact. For spin-one factors, a
submultiplicative row-sum norm proves the dimension-three trace bound without changing the matrices.

For the analytic layer, integrated Bessel inequalities give a basis-free Hilbert–Schmidt bound.
The rank-one expansion converges in operator norm, proving compactness without assuming a countable
basis. Kernel Hermitian symmetry proves self-adjointness through an explicitly integrable Fubini
calculation. These are uniform theorems, not numerical algorithms scaling with a seed size.

## Outcome

The actual complex-L² operator `K b` and its normalization `X a b = exp(-a*b) • K b` are constructed.
Both are Hilbert–Schmidt, compact, and self-adjoint. Their squared Hilbert–Schmidt bounds are
`9 * exp(6*b)` and `9 * exp((6 - 2*a)*b)`. Joint measurability, chronological order, null collisions,
empty-history value three, and symmetry are proved from the history and matrix definitions.

The operator trace and physical partition-function identification are separate results and are not
asserted here. No numerical spectral-gap certificate is produced by this prerequisite.
