# Finite seed redesign: candidate routes

All floating-point calculations are non-certifying. Exact scalar checks remain
conditional on the thermal enclosures, the trial inequality, and analytic transfer
lemmas; Julia is not the Lean trusted kernel.

1. Restrict existing source-1 data to even lengths through 10 or 8. Optimize
   nonnegative filters x^(2s) q(x)^2, including even q, and preserve the threefold
   N multiplicity using J=(Z+3P)/4 and N=(Z-P)/4. This avoids complex C traces.
   Cost is dominated by 3^n times the low-state count, with smaller polynomial
   degree reducing both table count and scalar proof complexity. Risk: the cap
   on J is too far above one, and its power grows at the trial length.
2. Compare the full I/O/N split, using all three twists, with the coarser J/N
   split and a twist-free bound. It may tighten caps enough to compensate for
   the extra Hamiltonian and complex arithmetic.
3. Raise temperature before extracting purities. Scan base beta, multiplier,
   moment order, and trial length. A longer trial improves spatial isolation but
   makes an imperfect upper cap more costly in the temporal estimate.
4. Change the energy shift only within a justified trial bound. Shifts rescale
   every transfer eigenvalue and moment; their benefit comes from slack in the
   variational inequality, not from changing a purity. Retain error rescaling.
5. Compare capped versus uncapped mass estimates, lower-degree filters, and
   relaxed absolute table errors. A cheap, slightly weaker seed is preferable
   to maximizing the gap if it removes the n=12 computations.

Chosen workflow: reuse exact stored thermal centers first, generate new small
ring spectra with Float64 magnetization-block ED/Lanczos only where necessary,
search scalar bounds, then replay 2–3 useful designs with exact rationals.
Report failures and conditional assumptions explicitly, along with work proxies
and the analytic lemmas each winning design needs. Do not infer a runtime for
an unimplemented Lean spectral checker from Julia timings.

## Outcome

Route 2 (I/O/N split with the C twist) won; route 1 (real twists only) works only with the
length-12 tables. Three designs replay exactly (README); the recommended one uses nine traces
at n = 6, 8, 10 with absolute error 3/10000 and gives a conditional gap > 0.005966152 for every
even L ≥ 60. Length 8 fails in the tested family even at zero error (best maximum defect about
0.529), and length 10 without the C twist fails as well (about 0.296). Route 4 (moving the
shift) only relabels the error budget. Details and the Pareto front are in METHOD.md.
