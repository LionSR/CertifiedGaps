# Square-lattice transverse-field Ising model: route analysis

The target is a uniform gap above a unique ground state on even square tori for
the Pauli Hamiltonian with ferromagnetic nearest-neighbor coupling and transverse field.
Both spatial periods must grow; a fixed-width ladder does not establish this statement.

## Candidate routes

1. **Product-state perturbation bounds.** At zero interaction the spatial purities are
   exactly one and the thermal purity is explicit. Finite-volume Weyl bounds give rigorous
   seeds at sufficiently large field, with negligible computational cost. Their error is
   proportional to inverse temperature times the number of bonds, so they are conservative.
   This provides a first complete certificate without assuming an infinite-volume gap.
2. **Verified finite-volume spectra.** The seed with periods two requires the spectra of
   the 2 by 2, 4 by 2, and 4 by 4 tori, at most sixteen spins. Existing core symmetry blocks
   and verified eigensolvers can enclose these spectra. The cost is exponential in seed
   area, but the inverse temperature only enters inexpensive partition sums. Period two
   includes repeated bonds, as required by the translation-invariant periodic sum.
3. **Strong-field connected-cluster expansion.** The bulk terms cancel in spatial
   partition-function ratios. Bounds on winding clusters could support larger seeds
   without diagonalizing larger tori. A useful implementation needs a proved convergence
   criterion and a quantitative remainder; an unverified series is insufficient.
4. **Classical spacetime transfer matrices.** A finite Trotter approximation has real
   symmetric spatial transfers in both directions. Its moment inequalities survive the
   Trotter limit. This supplies the structural input, but direct evaluation has a cost
   exponential in transverse width times Trotter depth.

## Chosen route

Prove the spatial scalar inequalities through finite Trotter transfers and their limit.
Use exact rational arithmetic for the three-component recurrence in the 2+1 dimensional
note. The verified sixteen-spin calculation proved inexpensive enough to proceed directly
with the spectral route; the product-perturbation estimate was not needed. No changes to
the manuscript or Lean formalization are required. The two-dimensional recurrence is
distinct from the core one-dimensional two-purity recurrence.

## Outcome

Complete verified spectra give certificates at $g/J=18,20,24$, with respective uniform
gap lower bounds $0.59437J$, $2.17363J$, and $4.80557J$ on every even square torus of side
at least four. All seeds have $n_0=m_0=2$ and $\beta_0g=4$. The largest spectral block
has dimension 2392. The full spectra take about 70 seconds per field on the development
machine with one BLAS thread.

The initial defects exceed $10^{-3}$. An exact three-component chain reaches that tail
region after five, four, and three steps, respectively. Bounding each earlier interpolation
window is necessary to retain coverage of all even square periods. The spatial-transfer
proof, seed enclosures, chain, and checker are documented in `METHOD.md` and `README.md`.

The numerical grid at $g/J=17$, with $\beta g=2:1/8:12$, found no closing seed. This is
a limitation of the tested finite seed and recurrence, not a gap exclusion. Further progress
toward smaller fields calls for larger seeds with controlled connected-cluster remainders
or sharper directional inequalities.

The row model exposed a generic cost in exact symmetry validation: dense multiplication
by signed permutation matrices. A separate minimal core change uses their sparse form,
preserving exact arithmetic and the same validation conditions. The new model uses the
shared core for all spectral reductions and eigenvalue enclosures.
