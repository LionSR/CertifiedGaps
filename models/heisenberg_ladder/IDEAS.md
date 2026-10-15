# Heisenberg ladder: route memo

Lengths count rungs (two spin-1/2 sites, d = 4 per cell). Survey values are hypotheses to
reproduce, never certified inputs. "Threshold" means the threshold reached by this
finite-size certificate method, not the physical gap-closing coupling.

## Candidate routes

1. **Closed forms.** At zero leg coupling Zₙ = (e^{3βJ/4} + 3e^{−βJ/4})ⁿ exactly. The
   interacting isotropic ladder is neither commuting nor free-fermionic (Jordan–Wigner leaves
   interactions), so this is only a reference point.
2. **Convergent expansions.** A strong-rung linked-cluster/polymer expansion could bound
   log Z₂ₙ − 2 log Zₙ exponentially in n at a cost set by the expansion order, not by n. At the
   couplings of interest a convergence majorant and a wrapping remainder would first have to
   be proved; none is available. High-temperature expansions do not reach the low-temperature
   seed. Not attempted.
3. **Transfer operator.** After the checkerboard rotation all cross-cell factors are real, so
   a Trotterised transfer kernel is positive, but its dimension grows exponentially with the
   number of slices and RP gives self-adjointness, not a contraction rate. Not attempted.
4. **Perturbation from independent rungs.** Weyl sandwiches log Zₙ around the rung product
   with an error extensive in n. Cost O(1), any n, but the bulk cancellation between the
   purities is lost. Tested analytically: **fails for every n and β when 0 < J⊥ ≤ 4**, even
   with the centred bond (PERTURBATION.md). The same norm bound is useful *locally in J⊥*: it
   continues a certified seed over |J⊥ − J*| ≤ δ.
5. **Smaller seeds.** P_τ(2n, β) = Pₛ(n, 2β) P_τ(n, β)² / Pₛ(n, β)² is exact, but bounding
   P_τ(2n) this way still needs Pₛ(n, 2β), i.e. Z₂ₙ(2β); no ring is saved without extra
   structure. Odd n₀ are not admissible, so ten-rung traces do not give an n₀ = 5 seed.
6. **Symmetry-resolved complete spectra (chosen).** Work in the unrotated SU(2) basis
   (unitarily equivalent to the RP model) with magnetisation, rung momentum, leg exchange and
   spin inversion; exact multiplicities from spin flip and complex conjugation. Dense cost is
   Σ O(block³), memory O(max block²). Verified enclosures (interval matrix, residual and
   orthogonality bounds, Weyl) for 2, 4, 8 rungs; floating-point scans through 10 rungs
   (largest block 8398) for orientation only.
7. **Variational sandwiches.** Rung-singlet or MPS trial energies and cluster lower bounds
   control E₀ but not the thermal tail, and a small energy-density error grows extensively
   in the purity ratio. Not attempted.

## Outcome

- Certified seeds at n₀ = 4 for J⊥ = 4, 7/2, 3, 5/2 and at the upper threshold endpoint, with
  coupling continuation (δ = 10⁻⁴) around 4, 7/2, 3.
- Every one of the 128 grid seeds n₀ ∈ {2, 4}, β₀ = k/4 is rigorously rejected (lower
  defects plus an exact recurrence witness) at J⊥ = 326133/131072; at J⊥ = 2 and 3/2 the
  rejection is extended to every β > 0 by monotone window bounds and a ground-energy tail.
- The floating-point scans through ten rungs find no even seed at J⊥ ≤ 2. A two-point fit of
  the spatial mismatch suggests n₀ ≈ 6 (J⊥ = 2) and 8 (J⊥ = 3/2), i.e. complete traces of
  24 and 32 spins, beyond dense ED. Reaching them needs route 2 or 3, not larger ED.

## Practical notes

- Running the six ten-rung scans concurrently with dense complex intermediates exceeded a
  17 GiB memory limit (one process OOM-killed). The numerical path therefore builds blocks
  sparsely and makes complex momentum blocks real with a sparse antiunitary change of basis.
  Run scans one or two at a time.
- The ten-rung spectra are stored compacted (sector table, log Z grid, lowest 32 levels);
  full level lists are regenerable with `scripts/scan.jl` and are not needed by any check.
