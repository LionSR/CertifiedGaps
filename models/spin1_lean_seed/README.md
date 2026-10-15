# Toward a Lean-checked seed for the spin-1 Heisenberg gap

The Lean proof of the spin-1 Heisenberg gap (`lean/CertifiedGaps/Heisenberg/Main.lean`) is
complete except for two seed inequalities. This directory collects four studies of how to prove
a seed inside the Lean kernel: with `decide +kernel` only, no `native_decide` and no new
axioms. **No seed and no gap is certified here.** Each subdirectory states what it checks
exactly and what is floating-point evidence.

The plan under study: Julia is an untrusted generator, and Lean checks every inequality.

- **Thermal traces.** Lower bounds come from trial vectors (Peierls). Upper bounds come from
  Rayleigh–Ritz plus eigenvalue lower bounds.
- **Eigenvalue lower bounds.** These need a *counting anchor* ρ ≤ λ_{k+1}. It comes from an
  operator comparison H ≥ H₀ that is checked locally.
- **Arithmetic.** Whole state vectors are packed into single big `Nat`s, because the kernel
  does those operations through GMP.

| study | question | main finding |
|---|---|---|
| [`seed_design/`](seed_design/) | which finite inputs does the seed need? | Three exact rational designs, each passing `ChainCertificate.Valid` and each conditional on the finite trace tables and the analytic seed lemmas. **Design C** uses even lengths 4–12, only the real twists 1 and P, and table error 10⁻³; it gives gap > 0.0390 for every even L ≥ 60. **Design A** uses lengths 6, 8, 10 with twists 1, P, C and error 3·10⁻⁴; it gives gap > 0.0060 for every even L ≥ 60. Length 8 alone fails, and so do the tested length-10 designs without C. |
| [`spectrum/`](spectrum/) | how many low states, how precise, how expensive? | At n = 12, β = 49/4 only 26 (twist 1) and 34 (twist P) states lie below the tail cutoff. 48-bit fixed-point Chebyshev filtering of degree 16 meets every residual target. Compressed Sᶻ indexing needs thousands of shift masks (3,816 at n = 12, Sᶻ = 0), so keep full 3ⁿ vectors, which need 2n masks. |
| [`anchor/`](anchor/) | can a cheap comparison anchor the eigenvalue count? | The local four-site PSD comparison works in exact arithmetic. It leaves too many comparison levels: 121 at n = 12, Sᶻ = 0, or 49 with a binned tail, against 8 physical levels. A one-shot anchor does not close. Float64 continuation needs 25–73 stages per sector, 658 in total. |
| [`lean_prototype/`](lean_prototype/) | does packed arithmetic fit the kernel? | Yes. One n = 12 matrix–vector product plus a dot product checks in 20–26 s (M5). Ritz/residual witnesses for H₆, H₈, H₁₀ check with the standard axioms only. Mathlib lacks Courant–Fischer, the Weyl inequalities and Cauchy interlacing. |

A fifth study, [`../lean_seed_small/`](../lean_seed_small/), carries the same pipeline furthest
on the cheapest model, the longitudinal Ising chain at g = 1/4, h = 1/20.

- **Kernel-checked arithmetic.** It checks the exact arithmetic for its 6- and 12-site rings
  with `decide +kernel`.
- **Proved bridge.** It proves the seed of `IsingLongitudinal` for that point from four
  explicit log Z enclosures.
- **Still open:** the min–max anchor, the scalar exp/trace bounds, and identification of the
  checked integer matrices with the Lean Hamiltonian.

## Where this leaves the plan

1. **Target design C.** It needs real arithmetic only, gives the strongest conditional bound,
   and its largest ring is 12.
2. **The bottleneck is the anchor.** A naive kernel continuation at n = 12 is estimated at tens
   of hours (`anchor/`), against a 15–60 min target. Cheaper anchors (better comparisons, fewer
   stages, Lehmann–Goerisch with a priori bounds) are the next research question.
3. **Formalization still needed.** The missing min–max/Weyl/Lehmann linear algebra, the
   twisted-sector and filter lemmas, and the bridge from packed base-3 matrices to the Lean
   Hamiltonian. The studies estimate 8–15k Lean lines for the whole seed route.

Each subdirectory follows the repository's directory format and has its own `check.jl` and
`test/runtests.jl`.
