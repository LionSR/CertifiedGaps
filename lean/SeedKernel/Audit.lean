import SeedKernel

/-! Axiom audit of the SeedKernel library. Every line must report a subset of `propext`,
`Classical.choice`, `Quot.sound` (or no axioms); `sorryAx` or `Lean.ofReduceBool` would mean a
`sorry` or a `native_decide`. `models/lean_seeds/check.jl` parses this output. -/

-- Shared packed lanes
#print axioms SeedKernel.radixPack_extract
#print axioms SeedKernel.pack_eq_radixPack
#print axioms SeedKernel.signed_floor_case
#print axioms SeedKernel.reverse_case
#print axioms SeedKernel.gram_case
#print axioms SeedKernel.dot_carry_counterexample
#print axioms SeedKernel.fast_reverse_case
#print axioms SeedKernel.gram_prepared_case

-- Spin-1 Heisenberg
#print axioms SeedKernel.Spin1.weyl_counting_one
#print axioms SeedKernel.Spin1.residual_anchor_not_enough
#print axioms SeedKernel.Spin1.h6_naive_anchor_separation_fails
#print axioms SeedKernel.Spin1.signed_h_cases
#print axioms SeedKernel.Spin1.filter_cases
#print axioms SeedKernel.Spin1.h6_ritz_witness
#print axioms SeedKernel.Spin1.h6_lane_bound
#print axioms SeedKernel.Spin1.h8_ritz_witness
#print axioms SeedKernel.Spin1.h8_lane_bound
#print axioms SeedKernel.Spin1.h10_ritz_witness
#print axioms SeedKernel.Spin1.h10_lane_bound

-- Longitudinal Ising
#print axioms SeedKernel.Ising.check6
#print axioms SeedKernel.Ising.check12
#print axioms SeedKernel.Ising.direct6
#print axioms SeedKernel.Ising.direct12_chunked
#print axioms SeedKernel.Ising.allChunks12_eq_range
#print axioms SeedKernel.Ising.direct12_canonical
#print axioms SeedKernel.Ising.coverage12
#print axioms SeedKernel.Ising.anchor6
#print axioms SeedKernel.Ising.anchor12
#print axioms SeedKernel.exists_eigenvalue_near
#print axioms SeedKernel.card_eigenvalues_lt_le
#print axioms SeedKernel.Ising.hamiltonian_mulVec_apply
#print axioms SeedKernel.Ising.residualSq_eq
#print axioms SeedKernel.Ising.lowSpectrum6
#print axioms SeedKernel.Ising.lowSpectrum12
#print axioms SeedKernel.Ising.seed_spatial
#print axioms SeedKernel.Ising.seed_thermal
#print axioms SeedKernel.Ising.first_point_seed
#print axioms SeedKernel.Ising.ising_gap
#print axioms SeedKernel.Ising.wrong_eigenvalue_rejected
#print axioms SeedKernel.Ising.wrong_vector_rejected

-- General seeds from low spectra (any Hermitian matrix, any RPRing and boundary)
#print axioms SeedKernel.le_re_quadratic_of_weights
#print axioms SeedKernel.exists_near_of_residual
#print axioms SeedKernel.lowSpectrum_of_near
#print axioms SeedKernel.thermal_ratio_ge
#print axioms SeedKernel.spatial_ratio_ge
#print axioms SeedKernel.trace_le
#print axioms SeedKernel.exp_le_expUpQ
#print axioms SeedKernel.spatialPurity_ge_of_check
#print axioms SeedKernel.thermalPurity_ge_of_check
#print axioms SeedKernel.openThermalPurity_ge_of_check
#print axioms SeedKernel.openTrace_le_of_check

-- Longitudinal Ising ring at any inverse temperature
#print axioms SeedKernel.Ising.ring_seeds_of_check
#print axioms SeedKernel.Ising.ring_gap_of_check

-- Free (open-boundary) longitudinal Ising chain
#print axioms SeedKernel.IsingOpen.count2
#print axioms SeedKernel.IsingOpen.weights2
#print axioms SeedKernel.IsingOpen.count12
#print axioms SeedKernel.IsingOpen.statSum12_1
#print axioms SeedKernel.IsingOpen.statSum12_2
#print axioms SeedKernel.IsingOpen.weightBlocks12_all
#print axioms SeedKernel.IsingOpen.openHamiltonian_mulVec_apply
#print axioms SeedKernel.IsingOpen.freeH_apply
#print axioms SeedKernel.IsingOpen.freeH_weights
#print axioms SeedKernel.IsingOpen.residualSq_eq
#print axioms SeedKernel.IsingOpen.weights_spec
#print axioms SeedKernel.IsingOpen.lowSpectrum2
#print axioms SeedKernel.IsingOpen.near8
#print axioms SeedKernel.IsingOpen.lowSpectrum12
#print axioms SeedKernel.IsingOpen.free_gap_of_check
#print axioms SeedKernel.IsingOpen.isingFreeCertificate_valid
#print axioms SeedKernel.IsingOpen.isingFreeCertificate3_valid
#print axioms SeedKernel.IsingOpen.isingFree_check
#print axioms SeedKernel.IsingOpen.isingFree_check3
#print axioms SeedKernel.IsingOpen.ising_free_gap
#print axioms SeedKernel.IsingOpen.ising_free_gap_long
