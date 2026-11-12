# Route memo: a small seed for spin-1

Leaderboard challenge `spin1-small-seed`: certify the spin-1 seed calculation on rings of at
most 12 sites to an error of $`10^{-3}`$, giving $`\Delta_L > 0.039`$ for all even $`L \ge 60`$,
small enough for Lean 4 to check.

## Candidate routes

1. **Dense exact diagonalization of the seed rings.** Rejected: the
   [`blbq_haldane`](../blbq_haldane/README.md) survey shows spin-1 rings reachable by dense ED
   ($`n \le 14`$) have purities $`\le 0.33`$, far below the admissible region
   $`\max(\varepsilon_{\mathrm{s}}, \varepsilon_\tau) < 0.12731`$. The seed scale is
   $`n_0 \approx 20\text{–}24`$ unshifted; the certified seeds live at $`n_0 = 36, 60, 2304`$
   via short-ring traces, not via diagonalizing the seed ring itself.
2. **Reuse the frozen spin-1 trace tables.** The spin-1 model
   ([`../spin1`](../spin1/README.md)) already computes, in exact integer/rational arithmetic,
   the 18 twisted traces on rings $`n = 4\text{–}12`$ at $`\beta = 49/4`$ (twists $`I, P`$),
   hash-pinned in `data/models/spin1/source1.json` with rational intervals of half-width
   $`3 \times 10^{-5}`$, audited against the upstream preprint's printed integers. These
   intervals certify the challenge's $`10^{-3}`$ tolerance with two orders of magnitude to
   spare. Cost: none — a hash audit plus scalar replay.
3. **Recompute the traces inside this model.** Rejected for now: the Poisson trace recursion is
   the spin-1 model's 37-minute `--full` job; duplicating it here adds no rigor because the
   frozen outputs are hash-audited. The from-scratch replay remains available at
   `code/models/spin1/check.jl --full`.
4. **Seed design on top of the tables.** The
   [`lean_seeds`](../lean_seeds/README.md) `seed_design` study already replayed three designs
   in exact rational arithmetic down to `ChainCertificate.Valid`. Design C (`real_even`) uses
   only even lengths $`n \in \{4,6,8,10,12\}`$ with the real twists $`I, P`$ at $`\beta = 49/4`$,
   lifts to $`\beta_0 = 49/2`$ (multiplier $`r = 2`$) at moment $`k = 30`$, and uses the stored
   length-120 integer MPS trial: defects $`p = 0.047914768`$, $`q = 0.112578548`$, and the
   9-step rational chain gives $`\Delta_L > 0.039036761`$ for even $`L \ge 60`$. Designs A/B use
   the complex $`C`$ twist and smaller gaps ($`\approx 0.006`$); C is the real-arithmetic,
   Lean-checkable one and dominates the published $`\log(20)/784`$ target.

## Chosen route

Route 2 + 4: hash-audit the frozen spin-1 trace intervals and trial, then replay Design C
entirely in `Rational{BigInt}` arithmetic (sector inversion $`J = (Z + 3P)/4`$, $`N = (Z-P)/4`$;
degree-4 squared filters with exact exterior rays; temperature lift; 9-step doubling chain on
the $`10^{-512}`$ rounding grid). No floating-point number accepts an inequality.

## Outcome

`certify.jl` writes `results/certificate.json`: seed $`(n_0, \beta_0) = (60, 49/2)`$ with exact
defects $`\varepsilon_{\mathrm{s}} = p`$, $`\varepsilon_\tau = q`$ (both below the diagonal
threshold $`0.12731`$), chain rate $`r = 130115651/50000000`$, and
$`\Delta_L > 39036761/10^9 = 0.039036761`$ for every even $`L \ge 60`$. `check.jl` (17 s) and
`check.jl --full` (18 s) replay everything and fail loudly on any mismatch; the test suite
(26 tests, 23 s) passes. The Lean 4 kernel check of the trace tables remains open (tracked by
`code/models/lean_seeds`).
