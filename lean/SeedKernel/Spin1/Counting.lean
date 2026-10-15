import Std
namespace SeedKernel.Spin1

/-- One-Ritz-value counting step in scaled integer coordinates. Weyl and the counted
    anchor are explicit premises. This is NOT an eigenvalue theorem for a Hamiltonian. -/
theorem weyl_counting_one (theta complement lambdaLow lambdaHigh rho eps : Int)
    (hWeylLow : min theta complement - eps ≤ lambdaLow)
    (hWeylHigh : lambdaHigh ≤ max theta complement + eps)
    (hAnchor : rho ≤ lambdaHigh)
    (hSeparation : theta + eps < rho) : theta - eps ≤ lambdaLow := by
  have h : theta ≤ complement := by omega
  simpa [Int.min_eq_left h] using hWeylLow

-- H=diag(-100,0), V=e₂ has zero residual. rho=0 bounds the second eigenvalue,
-- but theta+eps<rho fails. Residual and anchor alone do not give lambdaLow >= theta.
theorem residual_anchor_not_enough :
    min (-100 : Int) 0 = -100 ∧ max (-100 : Int) 0 = 0 ∧ ¬ (0 : Int) ≤ -100 := by
  decide +kernel

-- The naive cut anchor for H6 cannot isolate even the first measured Ritz value:
-- its second level is -10 (two open 3-site triplet ground spaces).
theorem h6_naive_anchor_separation_fails :
    ¬ ((-8495447 : Int) + 274967 < -10000000) := by decide +kernel

#print axioms weyl_counting_one
end SeedKernel.Spin1
