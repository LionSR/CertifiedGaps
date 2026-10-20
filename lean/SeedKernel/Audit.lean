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
#print axioms SeedKernel.Ising.exists_eigenvalue_near
#print axioms SeedKernel.Ising.card_eigenvalues_lt_le
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
#print axioms CertifiedGaps.SpinOneHeisenberg.trial_check_72
#print axioms CertifiedGaps.SpinOneHeisenberg.trial_check_120
#print axioms CertifiedGaps.SpinOneHeisenberg.gap_of_table
#print axioms CertifiedGaps.SpinOneHeisenberg.gap_log_twenty_of_table
#print axioms CertifiedGaps.SpinOneHeisenberg.seed60_of_table
#print axioms CertifiedGaps.SpinOneHeisenberg.gap_sixty_of_table
#print axioms CertifiedGaps.SpinOneHeisenberg.gap_sixty_chain_of_table
