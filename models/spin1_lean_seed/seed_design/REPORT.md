# Finite seed redesign for the Lean spin-1 gap proof

## Recommendation and trust boundary

**Use the nine source-2 traces at lengths 6, 8, 10 and twists `1, P, C`, allow absolute
error `3/10000`, retain the existing length-120 integer MPS, and lift `21/2` to `105/4`.**
Two degree-10 filters, no O-sector filter, and a length-30 intermediate moment give exact
rational seed bounds

- spatial defect at `(n0, beta) = (60, 105/4)`: `p < 0.124048525`;
- thermal defect at `(120, 105/4)`: `q < 0.125580476`;
- the rational chain implies **gap `> 0.005966152` for every even `L >= 60`**.

Both defects are below `0.12731`. The nine thermal dimensions sum to 199,017, versus
1,594,242 for the original eighteen inputs: **8.01 times smaller**. Including the measured
low-state counts gives a 9.93-fold reduction in the requested `sum 3^n K` proxy. Charging
the complex C computation two to four times a real trace still gives a 5.3–7.7-fold reduction.

These are **exact rational scalar checks conditional on physical finite inputs and analytic
lemmas**, not an unconditional Lean proof. Julia/Nemo is an untrusted generator. No Lean file
was changed, no `native_decide` was used, and no claim is made that a Julia assertion proves a
Lean theorem. The existing stored thermal enclosures were reused and hash-audited, not freshly
contracted. New ED/Lanczos data and all workload predictions are explicitly non-certifying.

The new route proves a stronger version of the final theorem via `ChainCheck.lean`; it is
**not a drop-in proof of the old two seed statements at exactly `(2304, 784)`**. The future
Lean change should route the final theorem through the new chain, then remove the obsolete
seed lemmas. The required spatial moment inequalities already exist in
`lean/CertifiedGaps/Heisenberg/Main.lean`; the twisted sector/filter input still needs a proof.

## Three exact rational replays

All three use the original shift `a = 700741/500000`. `b` is the temperature of the finite
thermal tables; `beta = r b` is the temperature of the resulting seed. Bounds below come from
`certificate.json`, with upward/downward rational rounding throughout.

| Design | Thermal lengths; twists | eps per trace | b | r; beta | Trial; n0 | p upper, displayed | q upper, displayed | Exact gap lower bound |
|---|---|---:|---:|---|---|---:|---:|---:|
| **A: small_simple, recommended** | 6,8,10; 1,P,C | 0.0003 | 10.5 | 5/2; 26.25 | 120; 60 | 0.124048525 | 0.125580476 | **0.005966152** |
| B: small_relaxed | 6,8,10; 1,P,C | 0.0004 | 10.5 | 25/8; 32.8125 | 144; 72 | 0.117346685 | 0.128925715 | **0.006309519** |
| C: real_even | 4,6,8,10,12; 1,P | 0.001 | 12.25 | 2; 24.5 | 120; 60 | 0.047914768 | 0.112578548 | **0.039036761** |

The covered lengths are **every even L at least n0**, with no claim about smaller lengths.
B's thermal defect exceeds the diagonal sufficient threshold; its full unequal-defect
`ChainCertificate.Valid` conditions pass. A and C satisfy both diagonal thresholds directly.
All three dominate the original `log(20)/784` target. This comparison was checked without
floating-point logarithms: for chain rate R and `784/beta = u/v`, the checker verifies
`R^u > 20^v`. For A this is `R^448 > 20^15`, where `R = 116954119/100000000`.

A deliberately uses `k=30`, hence spatial residual exponents 2 and 4. The numerical optimum
near `k=28` gains about 0.0007 in the bound but introduces fifteenth/seventh powers and roots.
B also uses `k=30`, with `h=12/5`; it trades more scalar proof work and a longer trial for a
small reduction in predicted low-state work. C uses `k=30`, `h=2`, and integer temperature
multiplier 2. These are proof-simplified choices, not claims of global numerical optimality.

### Filters and exact error propagation

The selected filters are `F(x)=x^s q(x)^2`. Each listed terminating decimal below denotes an
exact rational. The checker verifies the complete exterior rays by expanding `F(±(v+y))`
and checking every nonconstant coefficient is nonnegative, with the constant term strictly
above the summed moment bound.

| Design; sector | s | q(x) | Modulus cap v |
|---|---:|---|---:|
| A; I | 6 | `x^2 - 0.204792` | 1.0024667044 |
| A; N | 6 | `x^2 - 0.25255` | 0.8865462348 |
| B; I | 6 | `x^2 - 0.201511` | 1.0024860438 |
| B; N | 6 | `x^2 - 0.251916` | 0.8865829651 |
| C; J | 4 | `x^4 - 0.32529 x^2 + 0.008718` | 1.0018483311 |
| C; N | 4 | `x^4 - 0.386904 x^2 + 0.011397` | 0.8722877194 |

A and B use the uncapped O-sector mass. All selected filters use only even moments. Their
moment degrees are 10 for A/B and 12 for C. No odd-length thermal table is needed.

For the three-twist inversion,

```
I_n = (Z_n + 3 P_n + 8 C_n)/12
O_n = (Z_n + 3 P_n - 4 C_n)/12
N_n = (Z_n - P_n)/4.
```

Independent absolute trace errors eps imply sector errors `eps`, `2 eps/3`, `eps/2`,
respectively. For the two-twist inversion `J=(Z+3P)/4`, the errors are `eps` and `eps/2`.
The source code conservatively used eps for every sector; preserving these smaller errors
helps when loosening the thermal tables. The filter penalty is the exact coefficient
absolute sum times the corresponding sector error.

For A, the exact tenth-moment upper bounds are
`mI=1.042943416666...`, `mO=0.040230416666...`, `mN=0.40176425`.
Both capped sectors satisfy `v^10 < m < 2 v^10`. Thus only the two-piece mass bound is needed:

```
A_k = U^k + (mI-U^10)^(k/10)
B_k = mO^(k/10)
C_k = V^k + (mN-V^10)^(k/10)
M_k = A_k + 2 B_k + 3 C_k.
```

At k=30 and trial length 120 these moment exponents are integers. Raise the physical
Boltzmann factors by `r=5/2` to bound the new spatial moment by `D >= M_30^r` and both twisted
moments by `E >= (A_30+2B_30)^r`. Then
`tauI=(D+11E)/12`, `tauJ=(D+3E)/4`. The trial provides a fourth moment greater than one,
forcing a dominant I-sector weight. The residual convex bound is evaluated at the four
vertices of `0<=x<=tauI-t`, `0<=y`, `x+y<=tauJ-t`, with multiplicity divisors 2 and 3.
The physical residual ratio is bounded independently by `(M_120-1)^r`.

The replays use 12-decimal rational upper witnesses for D, E, Q and lower witnesses for the
dominant weight. Fractional powers are checked by integer-power comparisons, not floating
root evaluation. Nine rounded rational doubling steps use grid `10^-512`, then the existing
C=5 tail conditions. The three compact replay recipes together occupy **9,312 bytes**;
the expanded audit certificate occupies about 85 KiB. These sizes exclude the generic checker
and any future finite spectral certificate.

## Search coverage and numerical Pareto front

The search is finite and deliberately described precisely. It is not an impossibility proof
for arbitrary filters, moment inequalities, or trial tensors.

- Existing data: source 1 at `b=49/4`, source 2 at `b=21/2`; nmax 8, 10, 12 where available;
  Z-only, J/N and I/O/N decompositions; even-only and odd-inclusive filters; capped and
  uncapped masses; errors `3e-5,1e-4,3e-4,1e-3,3e-3,1e-2`.
- New thermal discovery: all three twists, lengths 4 through 10; b from 4 to 20 in steps of
  0.5, plus 12.25; errors `3e-5,1e-4,1e-3`; trial lengths
  `48,60,72,84,96,108,120,144,168`; even k from nmax through N/2; r from 1 to 6 in steps of 0.25.
- Filters: `x^s q(x)^2`, even s at least 4, degree at most nmax. Candidates come from regularized
  moment-matrix solves at seven endpoint guesses and five ridge values. Both even q and
  unrestricted q were considered. This is a search family, not all nonnegative polynomials.
- Source-2 refinement: errors through `0.00075,0.001`, independent removal of N and O caps,
  N through 192, and r in steps of 1/8. The final exact designs drop the O cap.
- Adversarial failures: b from 6 to 16 in steps of 0.25; N from 40 to 192 in steps of 4;
  r from 1 to 12 in steps of 0.25; zero-error optimistic cases; direct unlifted seeds too.
  Lift cases with `M_N>=2` are discarded: the resulting q bound is at least 3/4 for all r>=1,
  so this bound cannot enter the required bootstrap region.
- Shift scan: a displaced by `-0.02,-0.01,-0.005,0,+1e-7` from the original shift, retaining
  the rescaled physical and spatial trial lower bounds. The MPS energy was screened at
  every multiple of four from 24 through 192; only length 144 was freshly contracted exactly.

`pareto.csv` contains **39 nondominated rows from 94 successful grid representatives**, using
four objectives: minimize `sum 3^n K`, maximize eps, maximize the reported numerical gap,
and minimize the covered-length threshold. Representatives retain the best gap at each tested
trial length, rather than silently dropping shorter-length coverage. This is a Pareto front
of the recorded grid, at the canonical shift, not of all possible proofs. Every CSV row is
non-certifying, including rows built from the stored exact tables: their search arithmetic
and low-state cost estimates are floating point.

Selected numerical frontier points (the exact replays above sacrifice some gap for simpler
scalar proofs):

| nmax | eps | b -> beta | Traces | Trial; even L >= | Numerical gap estimate | Input status |
|---:|---:|---|---:|---|---:|---|
| 8 | even 0 | 6..16 -> up to 192 | up to 15 | 40..192 trials | none found | Float64/optimistic |
| 10 | 0.0005 | 10.5 -> 32.8125 | 9 | 144; 72 | 0.00332475 | stored tables |
| 10 | 0.0004 | 10.5 -> 32.8125 | 9 | 144; 72 | 0.00640352 | stored tables |
| 10 | 0.0003 | 10.5 -> 26.25 | 9 | 120; 60 | 0.00674278 | stored tables |
| 10 | 0.0001 | 10.5 -> 32.8125 | 9 | 144; 72 | 0.01191965 | stored tables |
| 10 | 0.0001 | 10.5 -> 31.5 | 21 | 144; 72 | 0.01784976 | new Float64 odd-twist inputs |
| 10 | 0.00003 | 10.5 -> 31.5 | 15 | 144; 72 | 0.01926514 | new Float64 odd-twist inputs |
| 12 | 0.003 | 12.25 -> 30.625 | 10 | 120; 60 | 0.02205770 | stored tables, even-only |
| 12 | 0.001 | 12.25 -> 29.09375 | 10 | 120; 60 | 0.04704958 | stored tables, even-only |
| 12 | 0.00003 | 12.25 -> 26.03125 | 18 | 120; 60 | 0.07210410 | stored tables, odd-inclusive |

The original experiment's optimized 0.072901771 bound uses a sharper filter and a bridge
covering L>=36. Our finite filter candidate search is not guaranteed to recover that optimum.
Do not replace that existing result with the slightly weaker last search row.

### Negative results that change the plan

**Length 8 does not look like a precision problem.** Even with zero error, the best observed
I-sector modulus cap is 1.006931714. Among 132,075 evaluated lift cases in the adversarial
zero-error grid, the best maximum defect was about 0.52877, at b=7.5, N=48, k=20, r=1.
With eps=3e-5 it was about 0.52902. There is no near-threshold length-8 winner in this family.
A successful nmax=8 route needs stronger moment information or a different analytic bound,
not just smaller numerical errors or another decimal in the current filters.

**Dropping C at length 10 fails in the tested family.** At the source-1 temperature, the best
J cap was about 1.00637. Over the finer zero-error beta grid it improved only to about
1.00449244; the smallest maximum defect among 478,980 evaluated lift cases was about 0.29558.
At eps=3e-5 the corresponding maximum defect was about 0.29664. This does not prove C traces
are logically necessary, but it kills the simple proposal to discard them while keeping these
short-moment bounds and trial scales.

**At canonical a, eps=0.001 is too loose for nmax=10 in the tested grid.** The finer scan
included 1,359,900 lift candidates. Its best maximum defect was about 0.1338344, with
`(b,N,k,r)=(10.75,148,30,3.25)`: both defects exceed the sufficient threshold, and the actual
unequal-defect recurrence fails as well. This is much closer to viable than the length-8 case.
The selected simpler source-2 designs work at eps=0.0003 and 0.0004.

**Keep the N cap, drop the O cap.** Removing the N cap found no successful source-2 seed in
the refinement grid even at eps=1e-5; removing only the O cap preserved useful seeds through
larger errors. Fully uncapped versions also failed in the existing-data search. The cap is
not dispensable merely because the high moments are small.

**A shift change does not create accuracy for free.** At a lowered by 0.02, eps=0.001 can
produce a numerical gap estimate about 0.0093664. But its effective original-normalization
error at n=10 is only `0.001 exp(-10.5*10*0.02) = 0.0001224564`.
The lower shift also raises the physical cutoff energy by `n*0.02`, requiring more low states.
The original length-120 trial permits increasing a by only about `1.0517e-7`, which rescales
a length-10 table by only about `1.1e-5` relatively. Keep the original shift unless the checker
has a separate normalization advantage; advertising the lowered-shift eps as cheaper would
misstate the workload.

## Twists, multiplicities, and analytic dependencies

**Yes: the threefold N multiplicity can be used without computing C-twisted traces.**
In `transfer.tex`, conjugation by the cyclic rotation C permutes the three nontrivial D2
character spaces, so their X restrictions have identical spectra. With J the trivial D2
space, `Z=J+3N` and `P=J-N`; identity and P traces therefore determine J and N.

One still needs the **C symmetry action/equivariance theorem** to prove those three spaces
are equivalent. Avoiding a C thermal computation does not remove that symmetry proof.
Computing C traces allows the finer split `J=I+2O`; complex conjugation pairs the two O
spaces. This is finite rotation-group representation theory, not a physical SU(2) multiplet
classification and not a momentum-block decomposition.

| Ingredient | A/B: nmax 10, I/O/N | C: nmax 12, J/N |
|---|---|---|
| C-twisted physical thermal trace | required by selected design | dropped |
| C action permuting N spaces | required | required |
| I/O splitting and twofold O multiplicity | required | dropped |
| Odd-length tables / odd spatial moments | dropped | dropped |
| Negative-eigenvalue sign isolation | dropped; only even moments | dropped; only even moments |
| O polynomial/cap | dropped | not present |
| Capped mass for the two controlled sectors | retained | retained |
| General floor-based capped-mass formula | can specialize to `h<m<2h` | same specialization |
| Dominant spatial weight and physical ground-weight arguments | both retained | both retained |
| Temperature lifting | retained | retained, multiplier 2 |
| Physical SU(2) or momentum/cyclotomic blocks | unnecessary | unnecessary |
| Periodic MPS contraction formula and variational principle | retained | retained |

A promising formalization simplification is to extend the **existing finite Trotter transfer
matrix** in `lean/CertifiedGaps/Ring/Transfer.lean`, rather than construct the paper's infinite
history-space operator. That file already proves `Z_(2N)` as a limit of moments of a finite PSD
matrix `T=F^T F`. All selected filters are even in x, so set `y=x^2`: A/B use
`y^3(y-c)^2`, C uses `y^2(y^2-cy+d)^2`. Prove the finite twisted trace identities and the
required rotation-sector equivalences for T, establish the finite scalar bounds with strict
margins, and pass to the Trotter limit. This could avoid the compact-operator construction,
countable spectral lists, and odd-moment identities entirely.

That extension is **not proved here**. Equivariance of the checkerboard transfer matrix,
the correct seam twist, multiplicities, and passage of the strict margins through the limit
must be established. Existing untwisted log-convexity alone does not supply these facts.

## Workload measurements and revised estimates

For each trace, predict the number K of physical states below

```
rho = log(2 * 3^n / eps)/b - n*a.
```

If a certified counting anchor later establishes this cutoff, the omitted-state tail is at
most eps/2, leaving half the error budget for retained-state enclosures. Here K is estimated
from **non-certifying** ED/Lanczos, not proved by the residuals. `summarize.jl` rejects a
truncated Ritz list that does not extend beyond rho. The estimates include physical
multiplicities and all Sz sectors; they do not assume SU(2) counting.

Predicted low-state totals for A, with the largest count in a single Sz block in parentheses:

| n | identity | P | C |
|---:|---:|---:|---:|
| 6 | 4 (2) | 6 (2) | 4 (2) |
| 8 | 9 (3) | 11 (4) | 9 (3) |
| 10 | 21 (7) | 18 (6) | 16 (6) |

These are far fewer states than a dense eigensolve, but their *completeness* is exactly the
unresolved counting-anchor obligation in the plan. Small Ritz residuals do not establish it.

| Design | Sum 3^n | Sum 3^n K | Sum d_m k_m | Sum d_m k_m^2 |
|---|---:|---:|---:|---:|
| Original 18 traces, eps=3e-5 | 1,594,242 | 34,237,485 | 4,553,325 | 27,198,539 |
| A | 199,017 | 3,448,170 | 489,700 | 2,393,682 |
| B | 199,017 | 3,211,974 | 462,640 | 2,339,562 |
| C | 1,195,722 | 21,977,649 | 2,879,341 | 14,620,013 |

The last two columns are coordinate-volume proxies for applying H to all retained vectors
and forming all Gram entries. They are not counts of packed GMP operations. They weight
negative/positive Sz copies by their physical multiplicity; an additional proven spin-flip
reduction could save work. Complex arithmetic factors are omitted in these raw columns.
With a factor-two C charge, A's state proxy is 4,454,919; with factor four it is 6,468,417.
The largest physical Sz block drops from 73,789 at n=12 to 8,953 at n=10.

The Float64 block eigensolves across lengths 4–12 took about **287 seconds of summed solve
time** on this host; matrix construction, compilation, and orchestration are additional.
At the two published temperatures, the new n<=10 traces differ from stored centers by at
most about `9.96e-7` (source 1) and `4.29e-7` (source 2), well within their original errors.
This is a useful independent control, not a proof of all Lanczos multiplicities or new traces.

The fresh **length-144 trial** used the 11 integer charge blocks of orders at most 162,
without SU(2) reduction. FLINT integer powers took about **5 seconds** (about 6 seconds with
result assembly); the exact energy per bond is approximately `-1.4014821051729465`.
The output contains the full U, V and positive integer margin `-500000 U - 700741 V`.
The three scalar replays take roughly 5–7 seconds including compilation; the two subsequent
warm replays take about 0.05 seconds each. Final tests passed **110 checks**, including rejection
of forged rates, invalid chain steps, inflated tails, and invalid filter caps.

### Lean effort and kernel time

The smaller n changes the **finite computation**, not the generic soundness proofs. A
reasonable revised planning decomposition, conditional on reusing the finite transfer API, is:

| Work item | Estimated new Lean LoC | Reason |
|---|---:|---|
| Packed integer action, Gram/residual identities, spectral/trace bounds and anchor interface | 2,000–3,200 | Generic checker still required; fewer n values do not remove its proof |
| C-twist exact arithmetic and I/O/N transfer-sector extension | 1,200–2,100 | Adds complex-pair arithmetic and finite rotation-sector proofs |
| Two simple filters, two-piece mass bound, isolation and temperature lifting | 600–1,000 | Only even moments; no O filter; simple A residual powers |
| Periodic MPS contraction and variational bridge | 600–1,000 | Independent of the thermal length reduction |
| Concrete witnesses and existing chain integration | 200–400 | Chain theorem is already present |
| **A total** | **4,600–7,700** | Estimate, not measured Lean development |

The real-only C route may save roughly 400–800 lines by avoiding C physical arithmetic and
the I/O split, while costing substantially more kernel computation. If a new infinite spatial
operator is required instead of extending the existing finite transfer, reserve another
1,000–2,500 lines; that is a risk allowance, not a theorem about implementation size.

The supplied four-million-bit multiplication benchmark does **not** calibrate packed H
applications, Gram extraction, or proof-term checking. Therefore an unconditional new total
kernel time cannot be inferred. If the original **15–60 minute thermal budget** is accurate
and scales with the complex-weighted state proxy, A projects to roughly **2–11 minutes for
that thermal component**. Add the unknown trial, analytic-proof, generation, and overhead
costs. For C the corresponding raw state-proxy projection is about **10–39 minutes**.
If the original 15–60 minutes included large fixed costs, even those projections are optimistic.
A measured kernel implementation must replace these conditional figures before promising a
total check time.

Do not retain the proposed 1.5–3 week schedule as an evidence-based commitment. At 250–400
reviewed Lean lines per working day the A estimate is roughly 3–6 engineer-weeks, assuming
the anchor succeeds and the finite transfer extension is straightforward. Neither assumption
was established by this task. The priority gate remains a sound, inexpensive low-spectrum
counting anchor at n=10; seed arithmetic is no longer the critical numerical obstacle.

## Deliverables, verification, and remaining risks

- `certificate.json`: three expanded rational scalar witnesses and chain certificates.
- `compact_witnesses.json`: 9.1 KiB of replay recipes for all three designs.
- `check.jl`, `test/runtests.jl`, `test_results.txt`: scalar replay and adversarial checks.
- `trial_scan.json`: fresh length-144 exact contraction plus explicitly non-certifying scans
  of other trial lengths. The exact source matrix recipe is checked against the pinned
  upstream manifest before contracting.
- `pareto.csv`, `summary.json`: the numerical grid frontier and workload estimates.
- `existing_search.json`, `new_search_8.json`, `new_search_10.json`, `refinement.json`,
  `diagnostics.json`, `shift_search.json`: successful and failed searches, with Julia scripts.
- `README.md`: reproducible commands. `Manifest.toml` pins package versions. Source and input
  hashes are recorded in outputs. The upstream paper was read at commit
  `adc7f1241b42e322a6451854ab7e4b4c146bf78a`, with downloaded files confined to this directory.

Remaining risks are concrete: the spectral counting anchor is unproved; a complex C trace
needs kernel-checked arithmetic; the finite twisted transfer/sector extension is missing;
Peierls/Ritz and MPS formulas must be tied to the Lean Hamiltonian and its normalization;
and the packed kernel's speed has not been measured on these witnesses. An exact scalar
recheck does not close any of those obligations. The nmax=8 and no-C nmax=10 failures are
negative results for the tested family, not universal lower bounds.

All authored and generated files are confined to `models/spin1_lean_seed/seed_design/`.
Neither `lean/` nor `models/spin1/results/` was modified.
