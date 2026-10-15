# Toward a Lean-checked seed: methods and findings

One section per study, each the full method document of the study as it was written (headings
demoted one level, paths updated to this directory's layout: `results/<study>/`,
`scripts/<study>/`, `src/<study>.jl`, and `lean/SeedKernel/` for the Lean sources). README.md
summarises the results; IDEAS.md has the route memos.

1. `seed_design`: finite seed redesign for the spin-1 gap proof.
2. `spectrum`: seed spectra, filtering precision and kernel costs.
3. `anchor`: comparison anchors, continuation and full tables.
4. `spin1_kernel`: packed-`Nat` prototype, Mathlib audit, mathematical review and budget.
5. `ising`: the two-level route for the longitudinal Ising seed and the 43-target inventory.

## 1. seed_design: finite seed redesign for the Lean spin-1 gap proof

### Recommendation and trust boundary

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

### Three exact rational replays

All three use the original shift `a = 700741/500000`. `b` is the temperature of the finite
thermal tables; `beta = r b` is the temperature of the resulting seed. Bounds below come from
`results/seed_design/certificate.json`, with upward/downward rational rounding throughout.

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

#### Filters and exact error propagation

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

### Search coverage and numerical Pareto front

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

`results/seed_design/pareto.csv` contains **39 nondominated rows from 94 successful grid
representatives**, using four objectives: minimize `sum 3^n K`, maximize eps, maximize the
reported numerical gap, and minimize the covered-length threshold. Representatives retain the best gap at each tested trial length,
rather than silently dropping shorter-length coverage. This is a Pareto front
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

#### Negative results that change the plan

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

### Twists, multiplicities, and analytic dependencies

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

### Workload measurements and revised estimates

For each trace, predict the number K of physical states below

```
rho = log(2 * 3^n / eps)/b - n*a.
```

If a certified counting anchor later establishes this cutoff, the omitted-state tail is at
most eps/2, leaving half the error budget for retained-state enclosures. Here K is estimated
from **non-certifying** ED/Lanczos, not proved by the residuals. `scripts/seed_design/summarize.jl` rejects a
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
warm replays take about 0.05 seconds each. Final tests passed **110 checks** (46 exact, plus
64 Float64 spectrum controls that run only when `cache/seed_design/` holds the n = 4–10 spectra),
including rejection of forged rates, invalid chain steps, inflated tails, and invalid filter caps.

#### Lean effort and kernel time

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

### Remaining risks

The open obligations are concrete: the spectral counting anchor is unproved; a complex C trace
needs kernel-checked arithmetic; the finite twisted transfer/sector extension is missing;
Peierls/Ritz and MPS formulas must be tied to the Lean Hamiltonian and its normalization;
and the packed kernel's speed has not been measured on these witnesses. An exact scalar
recheck does not close any of those obligations. The nmax=8 and no-C nmax=10 failures are
negative results for the tested family, not universal lower bounds.

The upstream paper was read at commit `adc7f1241b42e322a6451854ab7e4b4c146bf78a`
(`scripts/fetch_paper.jl`); downloads stay in the ignored `upstream/` directory.

## 2. spectrum: spin-1 seed spectra, filtering precision and kernel costs

Long-form derivations and full tables behind README.md (originally `REPORT.md`, generated by
the since-removed `report.jl`). File names refer to `results/spectrum/`; `eigenvalues.csv` and
`cost_sectors.csv` are regenerated by `scripts/spectrum/analyze.jl` and `scripts/spectrum/summarize.jl`.

**All spectral, filtering, anchor, and runtime projections below are non-certifying.**
Only the standalone Lean arithmetic benchmarks were checked by the Lean kernel.
No seed proof, native_decide, sorry, or extra axiom was introduced. No file under
`lean/` or `models/spin1/results/` was changed.

### Findings that change the plan

The number of low states is manageable, and 64-bit lanes with short filtering
cycles can meet the numerical residual targets. The proposed combination of a
compressed Sz vector and approximately 50 masked shifts is unsupported: direct
compressed indexing needs thousands of displacement masks. The spectral anchor
is also unresolved; the elementary two-half-chain comparison is much too weak.
These are larger risks than obtaining sufficiently accurate Ritz vectors.

At n=12 and beta=12.25, the requested dimension-tail criterion needs a cutoff
about **2.20335 above E0**, retaining **26 states for I and 34 for P** including
magnetization multiplicities. At beta=10.5 the same n=12 criterion needs about
2.57 and 58/60 states. Source 2 stops at n=10: its three twists need cutoffs near
2.36 and retain 22, 30, and 26 states for I, P, and C.

For n=12, I, Sz=0, a degree-32 Q48 recurrence plateaued at residual **2.99e-6**
after 40 cycles. Q112 reached **2.34e-9** in four cycles. Reducing the polynomial
to degree 16 recovered Q48: **3.29e-8** in eight cycles. More precision is not
automatically the cheapest repair; frequent orthogonalization controls the loss
of small filtered components.

Recommended next experiment: implement one exact, replayable Q48 degree-16
generator and Gram checker for n=10, including its rational mixing recipe, then
demonstrate a useful certified anchor. Prefer source 2 for the first prototype.
Defer a full Lehmann–Goerisch development until the anchor is known to work.

### Scope, conventions, and validation

Computed 21 ring/twist pairs: n=4..12 with I/P, and n=6,8,10 with C. The 357
signed-sector rows include dimensions and low eigenvalues. The 189 nonnegative
sectors were actually diagonalized; negative sectors are marked as inferred.
Spin inversion maps I/P sectors m to -m; for C it reverses the twist, and complex
conjugation restores spectral equality. This requires only an elementary basis
map, not SU(2) multiplet theory.

The Hamiltonian has diagonal sum s_j*s_(j+1), off-diagonal coefficients 1, and
closing-bond coefficients 1, -1, or exp(±2pi*i/3). Bond eigenvalues are -2,-1,1,
so n is an analytic upper bound on the spectrum. The shift is a=700741/500000.

“Exact diagonalization” here means **full dense Float64 diagonalization**, not
exact arithmetic. All sectors through n=8 were fully diagonalized. Above n=8,
ARPACK requested 32 eigenpairs, expanding to 64 when needed, with guard levels
beyond E0+3.15. Explicit residuals were below 1e-8 (typically 1e-13). The
full-spectrum trace and second moment, positive/negative sector equality, n=4
ground energy -6, and dense/Lanczos multiplicities were checked. The independent
Lanczos comparison through n=8 passed 143 checks. Integer floor arithmetic passed
1,009 tests, including comparison against BigInt. Final data consistency passed
307 checks, including all 27 stored reference trace centers to 3e-5.

Residuals and repeated successful numerical comparisons **do not prove that
Lanczos omitted no eigenvalues**. Every tail and anchor count here remains
conditional on the observed spectral ordering. Existing trace certificates were
read only as cross-checks; no certified claim is imported from them.

The upstream `references/periodic-gap-source/` directory was absent in this
checkout. Conventions and numerical centers were checked against the existing
Julia implementation and result files. Julia 1.12.5, Arpack 0.5.4, and the full
package manifest are recorded locally. The repository started at commit
463c4a5b1e345c762133d4ee71bbe79943a7db47.

### Low-spectrum counts and tail cutoffs

Each k counts all Sz sectors with multiplicity. Values within 1e-8 of an exact
fixed cutoff are treated as boundary values and excluded; this matters at n=4.
Let Q=sum exp(-beta*(lambda-E0)). The reported c* is the numerical infimum for
`(3^n-k(c))*exp(-beta*c)/Q < 1e-6`. Use a strict positive margin above c* in an
actual certificate. The energy shift cancels from this relative criterion.

| n | g | E0 | k(1) | k(1.5) | k(2) | k(2.5) | k(3) | c* / k*, beta=12.25 | c* / k*, beta=10.5 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 4 | I | -6 | 1 | 4 | 4 | 4 | 4 | 1.48239 / 4 | 1.72945 / 4 |
| 4 | P | -5.38164 | 2 | 2 | 7 | 11 | 13 | 1.48239 / 2 | 1.72159 / 7 |
| 5 | I | -6.53113 | 7 | 10 | 10 | 10 | 28 | 1.57278 / 10 | 1.83489 / 10 |
| 5 | P | -7.09592 | 3 | 3 | 5 | 11 | 17 | 1.5751 / 3 | 1.8366 / 5 |
| 6 | I | -8.61742 | 4 | 4 | 4 | 15 | 21 | 1.66541 / 4 | 1.94287 / 4 |
| 6 | P | -8.30592 | 2 | 7 | 11 | 14 | 26 | 1.66415 / 7 | 1.94 / 11 |
| 6 | C | -8.40229 | 4 | 4 | 11 | 16 | 22 | 1.66518 / 4 | 1.94165 / 9 |
| 7 | I | -9.5722 | 7 | 10 | 10 | 26 | 44 | 1.75519 / 10 | 2.04766 / 10 |
| 7 | P | -9.86337 | 3 | 5 | 11 | 17 | 28 | 1.75512 / 5 | 2.0469 / 11 |
| 8 | I | -11.337 | 4 | 4 | 15 | 21 | 32 | 1.84498 / 9 | 2.15194 / 21 |
| 8 | P | -11.1562 | 6 | 7 | 14 | 24 | 50 | 1.84438 / 13 | 2.15089 / 18 |
| 8 | C | -11.2093 | 4 | 7 | 16 | 22 | 45 | 1.84475 / 11 | 2.15159 / 16 |
| 9 | I | -12.48 | 7 | 10 | 20 | 34 | 63 | 1.93482 / 20 | 2.25712 / 28 |
| 9 | P | -12.6468 | 5 | 5 | 17 | 26 | 49 | 1.9344 / 17 | 2.25608 / 20 |
| 10 | I | -14.0941 | 4 | 9 | 21 | 26 | 67 | 2.0242 / 21 | 2.36088 / 22 |
| 10 | P | -13.9821 | 6 | 11 | 20 | 36 | 74 | 2.02395 / 20 | 2.36048 / 30 |
| 10 | C | -14.0138 | 4 | 9 | 20 | 35 | 70 | 2.02412 / 20 | 2.36076 / 26 |
| 11 | I | -15.3367 | 7 | 10 | 28 | 49 | 87 | 2.1142 / 34 | 2.46631 / 49 |
| 11 | P | -15.438 | 5 | 15 | 22 | 47 | 90 | 2.11366 / 22 | 2.46502 / 43 |
| 12 | I | -16.8696 | 4 | 15 | 22 | 48 | 91 | 2.20334 / 26 | 2.56956 / 58 |
| 12 | P | -16.7976 | 6 | 14 | 24 | 56 | 99 | 2.20335 / 34 | 2.56974 / 60 |

At n=12 the c=2 dimension-tail bound is approximately 1.21e-5 relative at
beta=12.25 and 3.96e-4 at beta=10.5. Thus “about E0+2” understates the needed
cutoff at the largest length. `sectors.csv` gives each individual sector's counts;
`eigenvalues.csv` (regenerated by `analyze.jl`) gives every computed eigenvalue. Bounds on the omitted terms
using the last computed guard eigenvalues are in `traces.csv`; these too are
numerical, conditional on complete counting.

### Absolute trace tolerance and eigenvalue accuracy

Write Zhat=Tr exp(-beta*(H+a*n)). A central numerical spectrum with trace Ztilde
and eigenvalue intervals ±epsilon gives width `2*Ztilde*sinh(beta*epsilon)`.
Using the observed Zhat as the planning approximation to Ztilde gives:

- A full-width budget delta permits epsilon=asinh(delta/(2*Zhat))/beta.
- Absolute error (radius) delta permits epsilon=log(1+delta/Zhat)/beta.
- The filtering experiment uses the stricter epsilon=asinh(delta/(4*Zhat))/beta,
  reserving half of a **full-width** budget for eigenvalue uncertainty. It also
  requires the upper tail below delta/4. This is stricter than interpreting the
  user's 3e-5 solely as an absolute radius.

The following ranges are for absolute-radius delta and cover all entries of
each source table. Every entry, sector, beta, and tolerance is in `precision.csv`.

| Source | Absolute trace error | Smallest epsilon | Largest epsilon |
| --- | --- | --- | --- |
| 1 | 3e-05 | 1.96098e-08 | 0.000832942 |
| 1 | 0.0001 | 6.5366e-08 | 0.00274404 |
| 1 | 0.001 | 6.53658e-07 | 0.0240046 |
| 2 | 3e-05 | 3.19401e-07 | 8.2398e-06 |
| 2 | 0.0001 | 1.06467e-06 | 2.74632e-05 |
| 2 | 0.001 | 1.06461e-05 | 0.000274277 |

The tightest source-1 entry is n=4/I: Zhat=124.88538, so absolute error 3e-5
requires about 1.96e-8 uniform eigenvalue accuracy; our half-width allocation
uses 4.90e-9. At n=12/I the absolute-radius allowance is about 1.29e-6. Precision
does not become uniformly harder with n. Odd rings can have very small shifted
traces and much looser absolute requirements.

A 1e-6 *relative* tail is not always small enough in absolute terms: n=4/I would
allow about 1.25e-4. The code raises its working cutoff from 1.48239 to 1.71198
to meet the reserved absolute-tail budget. The retained state count is unchanged.
Rounding of exponentials and final rational interval arithmetic still need their
own certified error allocation; no exponential was computed inside Lean here.

### Filter experiments and residual targets

The block starts from distinct computational-basis product states chosen
deterministically from a broad low-diagonal pool. Its width is p=k+8, capped by
the sector dimension. An empirically selected passband endpoint lies between
the p-th and (p+1)-th eigenvalues. A normalized Chebyshev recurrence suppresses
the unwanted interval [alpha,n]. After each degree-16/32/64 cycle, QR and
Rayleigh–Ritz are performed. Every scalar-vector H product is counted, including
the p products for Rayleigh–Ritz after each cycle.

The scalar recurrence coefficients are rounded to 40 fractional bits. “64-bit”
means Q48 with 16 sign/headroom bits; “128-bit” means Q112. The recurrence uses
exact signed Int128 arithmetic and floor divisions (with a split multiplication
primitive that avoids overflow). C uses two integer Eisenstein coordinates.
Native overflow headroom was checked in the measured runs. **QR, eigensolves,
input/output conversion at QR boundaries, and residual evaluation remain
Float64.** Q112 recurrence arithmetic does not make the whole experiment a
112-bit spectral computation.

For the linear proxy the block operator residual must be <=epsilon. For the
quadratic proxy the target is sqrt(epsilon*g/2), where g=rho-lambda_k is the
empirical anchor margin. This is a conservative scaling model for a quadratic
enclosure, **not an implemented Lehmann–Goerisch theorem**. Both require a
certified count/anchor and correct eigenvalue association before they can bound Z.

All active sectors were measured at Float64 degrees 16/32/64 and Q48/Q112 degree
32. Additional Q48/Q112 degree-16/64 sweeps cover n=10 and n=12, including C at
n=10. This is not a full precision-by-degree grid at every n. The adaptive rows
below select the successful measured degree with the fewest H products in each
sector; `filter_runs.csv` records the chosen degree and `filter_curves.csv` records
every cycle, residual, and Ritz error.

| Source | Lane bits | Target model | Successful sectors | Total Hv | Scalar Gram products | Mixing recipe estimate, KiB |
| --- | --- | --- | --- | --- | --- | --- |
| 1 | 64 | linear | 44/44 | 34388 | 1977 | 269.812 |
| 1 | 64 | quadratic | 44/44 | 26562 | 1977 | 229.469 |
| 1 | 128 | linear | 44/44 | 33879 | 1977 | 210.359 |
| 1 | 128 | quadratic | 44/44 | 26562 | 1977 | 229.469 |
| 2 | 64 | linear | 24/24 | 17800 | 2100 | 183.312 |
| 2 | 64 | quadratic | 24/24 | 12486 | 2100 | 121.734 |
| 2 | 128 | linear | 24/24 | 17770 | 2100 | 172.25 |
| 2 | 128 | quadratic | 24/24 | 12439 | 2100 | 106.703 |

There are 44 active nonnegative sectors for source 1 and 24 for source 2. Sectors
with k=0 need no trial vectors but **still require a rigorous exclusion/anchor**.
Negative-sector equality is reused for cost counts. Without that elementary
symmetry proof, repeat the nonzero-m sectors.

Uniform degree-32 Q48 failed the linear target in two sectors in each source.
All those failures were repaired by measured shorter cycles. Q112 degree-32
passed every sector. Degree-64 Q48 can be much worse: n=12/P/Sz=0 stalled around
0.0666. Increasing degree without increasing precision can destroy the filtered
block's usable rank. The quadratic target reduces adaptive Q48 Hv counts by
about 23% for source 1 and 30% for source 2, while final Gram costs are unchanged.

**Recipe gap:** the prototype feeds numerical QR output back into the integer
recurrence; it does not replay QR through stored dyadic mixing matrices. The
recipe-size column budgets 16 bytes per real mixing coefficient, including both
coordinates for C, and counts p² coefficients per cycle. These are estimates,
not exported validated recipes. Integer replay, coefficient magnitude bounds,
and a stable exact Gram check must be demonstrated before accepting the kernel
cost forecast. Initial product seeds and polynomial coefficients alone do not
describe these measured multi-cycle runs.

### Packing and operation counts

For sorted compressed Sz basis indices, the number of distinct displacement
classes depends strongly on the sector. The following counts include both
off-diagonal directions per bond. Diagonal contributions are additional.

| n | Sz | Sector dimension | Full lanes | Compressed masks | Full masks |
| --- | --- | --- | --- | --- | --- |
| 8 | 0 | 1107 | 6561 | 126 | 16 |
| 8 | 1 | 1016 | 6561 | 268 | 16 |
| 8 | 2 | 784 | 6561 | 320 | 16 |
| 10 | 0 | 8953 | 59049 | 624 | 20 |
| 10 | 1 | 8350 | 59049 | 1514 | 20 |
| 10 | 2 | 6765 | 59049 | 1970 | 20 |
| 12 | 0 | 73789 | 531441 | 3816 | 24 |
| 12 | 1 | 69576 | 531441 | 9884 | 24 |
| 12 | 2 | 58278 | 531441 | 13620 | 24 |

At n=12/Sz=0 a Q48 compressed vector has 4,722,496 bits; a full ternary vector
has 34,012,224 bits. Q112 doubles those widths. Full ternary indexing has 24
off-diagonal shifts, while the compressed Sz=0 and Sz=2 bases need 3,816 and
13,620 direct displacement masks. A structured routing algorithm could improve
compressed indexing; none was implemented here. The scan disproves the direct
masked-shift assumption, not every possible compressed layout.

Materializing just the n=12/Sz=0 compressed displacement masks at 64 bits per
lane would take about 2.10 GiB. The 24 full-basis masks take about 97.3 MiB.
Streaming or formula-based mask generation can reduce stored memory, but its
cost and proof must be included; it is not part of a kilobyte-sized certificate.

`cost_sectors.csv` (regenerated by `summarize.jl`) and `cost_totals.csv` provide four models: full ternary;
full Hv plus one final compaction for Gram products; direct compressed shifts;
and an explicitly hypothetical compressed layout with full-basis shift counts.
The hybrid compaction is a costed proposal using a fixed permutation network,
not implemented packing code. Its mask construction and proof remain open.

Counting rules: three primitive passes per off-diagonal displacement (mask,
shift, add), binary diagonal masks, 11 recurrence passes per vector product, and
four passes per mixing coefficient. C accounts for two coordinates and closing
bond mixing. Gram checks use three k-by-k symmetric Gram matrices, with a
conservative four-real-product complex dot product. Packed dot-product lanes
have width 2w+ceil(log2(d))+4, preventing cross-lane accumulation in the model.
Signed encoding, carry handling, and extraction still need a verified design.

The “bit-pass volume” is sum(width * primitive passes); it is a traffic proxy,
not a literal CPU instruction count. General multiplication is counted separately
using M(B)=B*log2(B)*log2(log2(B)). Both Gram and combined model bit-operation
counts are in the CSV. Small-iteration overhead is charged at 30,000/s. These
are transparent forecasts, not measurements of an implemented packed checker.

| Source | Layout, Q48 linear | Hv | Bit-pass volume | Gram 4M-bit equivalents | Combined model bit operations |
| --- | --- | --- | --- | --- | --- |
| 1 | full_ternary | 34388 | 4.14673e+13 | 20146.1 | 4.93407e+13 |
| 1 | full_Hv_compact_Gram | 34388 | 4.29641e+13 | 2361.34 | 4.3887e+13 |
| 1 | compressed_direct | 34388 | 9.57938e+14 | 2361.34 | 9.58861e+14 |
| 2 | full_ternary | 17800 | 4.63732e+12 | 3031.33 | 5.82201e+12 |
| 2 | full_Hv_compact_Gram | 17800 | 4.8151e+12 | 375.178 | 4.96173e+12 |
| 2 | compressed_direct | 17800 | 2.432e+13 | 375.178 | 2.44667e+13 |

### Lean timings and revised runtime budget

Installed the pinned Lean 4.35.0-rc3 from its official release after elan's
download endpoint failed. The standalone four-million-bit product/division
identity took 0.91 s including startup. Bulk kernel-checked loops amortize
startup more reliably: 100 four-million-bit multiplications took 2.26–3.08 s
including startup, approximately 0.016–0.024 s per multiplication after subtracting
the warm baseline. Thirty 20M-bit products took 4.14–6.38 s, and ten 68M-bit
products 4.98–5.17 s. The operands are dense patterned integers, not powers of two.

Bulk mask/shift/XOR loops give approximately 5.6e-11 to 9.4e-11 seconds per bit
pass. The 30k-to-100k small-loop timing difference corresponds to about 31k
iterations/s, consistent with the supplied 30k/s premise. These are elapsed
times on this machine, with one other Julia process active during some runs;
they are not portable performance guarantees. `kernel_bulk.csv`,
`kernel_timings.csv`, and the Lean-generating Julia scripts preserve the evidence.
Printed axioms were either empty or `[propext]`.

The forecast uses 5e-11..1e-10 seconds per bit pass and 0.015..0.035 s per
4M-equivalent multiplication, plus the small-loop charge. The following estimates
cover **vector generation, modeled mixing, and final Gram arithmetic only**.
They exclude anchor construction/recursion, scalar enclosure proofs, analytic
transfer lemmas, MPS trials, mask generation, and proof-term memory overhead.

| Source | Bits | Layout | Measured-rate model, min | Using 1 s / 4M product, min |
| --- | --- | --- | --- | --- |
| 1 | 64 | full_ternary | 41.2553–82.5267 | 406.543 |
| 1 | 64 | full_Hv_compact_Gram | 38.1192–74.7097 | 112.688 |
| 1 | 64 | compressed_direct | 959.392–1758.46 | 1796.44 |
| 1 | 128 | full_ternary | 77.1525–155.949 | 787.167 |
| 1 | 128 | full_Hv_compact_Gram | 71.0554–140.807 | 215.245 |
| 1 | 128 | compressed_direct | 1712.66–3268.26 | 3342.7 |
| 2 | 64 | full_ternary | 5.72776–10.6026 | 59.3566 |
| 2 | 64 | full_Hv_compact_Gram | 5.25516–9.3928 | 15.4269 |
| 2 | 64 | compressed_direct | 49.4751–69.8669 | 75.901 |
| 2 | 128 | full_ternary | 10.2843–19.9665 | 116.346 |
| 2 | 128 | full_Hv_compact_Gram | 9.31262–17.5426 | 29.5619 |
| 2 | 128 | compressed_direct | 69.646–110.311 | 122.33 |

For planning, allow **1–4 hours for source-1 thermal vector/Gram checking** and
**10–30 minutes for source 2**, conditional on a working Q48 replay and efficient
full-index packing. These are engineering allowances over the measured-rate
models, not upper bounds for the complete proof. A direct compressed-shift
checker is substantially worse. The supplied <1 s multiplication fact alone
supports only the much larger final column; it cannot justify 15–60 minutes.

The smaller source-2 table reduces Q48 full-layout bit-pass work by roughly a
factor of nine even though its lower temperature requires more states in some
sectors. This supports n_max reduction as a major lever. It does not establish
that arbitrary n_max=8 filter or seed redesigns preserve the desired gap.

### Anchor probe: the unresolved bottleneck

The exact local inequality h>=-2I gives a simple comparison obtained by cutting
two bonds: H >= H_open(floor(n/2)) + H_open(ceil(n/2)) - 4I. Full dense numerical
half-chain spectra through length 6 were convolved by Sz to test this baseline.

At n=12/I/Sz=0, the required anchor is rho=-14.6662193, with k=10. The naive
comparison has **242** eigenvalues below rho; its 11th eigenvalue is about
**-17.0810228**, a **2.4148** energy deficit. Across all sectors it has **968**
levels below the anchor, versus the ring's retained 26. It also fails to exclude
the k=0 Sz=3 and Sz=4 sectors. At n=10/C/Sz=0 (source 2), the corresponding
numbers are k=8 versus 165 levels and an anchor deficit about 2.5824.

`anchor_probe.csv` gives every source/sector result. Tuned boundary terms or a
homotopy comparison may improve this; this experiment did not test them. It is
incorrect to turn a small Lanczos residual into a certified upper trace without
closing this counting gap. The naive anchor fails; the tuned proposal remains
an open problem rather than a confirmed implementation detail.

### Revised implementation estimate and open risks

The original 4.5–7.5k Lean-line estimate cannot be validated from numerical
experiments. A working engineering budget is **8–15k Lean lines** for the whole
proposed route, with these additions to the original baseline:

| Additional work | Rough Lean lines | Reason |
| --- | ---: | --- |
| Full-index packing, signed lanes, optional compaction | 1–2k | Compressed fixed-shift premise fails |
| Replayable mixers and exact Gram conditioning/counting | 1–2k | Numerical QR is not a small kernel recipe yet |
| Recursive/tuned anchor or homotopy infrastructure | 1.5–3k | Elementary split comparison misses the count badly |

These are judgment ranges, not measured code-generation rates, and overlap may
reduce them. Deferring quadratic eigenvalue theory could avoid another 1–2k
lines; the measured Hv savings alone do not clearly justify it. The original
1.5–3 week completion claim is unsupported while the anchor is unresolved.
Budget 3–6 weeks after a successful anchor/replay prototype, with an unbounded
research risk if that prototype fails; this is not a delivery commitment.

Outstanding gates are: a kernel-checkable count anchor in every sector; exact
replay of the measured orthogonalization schedule; Gram conditioning and signed
packed arithmetic bounds; actual full-index Hv/Gram benchmarks at target sizes;
and a complete absolute-error budget including rational exponential enclosures.
The scalar transfer/MPS portions of the original plan were not investigated in
this task. There is no unconditional gap proof in this deliverable.

### Artifacts and reproduction

See README.md, "Reproduce" and "Layout".

## 3. anchor: comparison anchors for the spin-one seed, with full tables

This is the long companion to [README.md](README.md): the comparison construction, what the
counting anchors do and do not prove, the binned tail, the continuation that was tested, the
rejected routes, and the Lean cost estimate. The tables after the marker at the end are
regenerated by `scripts/anchor/report_tables.jl`. Everything labelled NON-CERTIFYING is Float64.

### The comparison and its exact check

Use the integer spin-one bond `h = Sz⊗Sz + (S+⊗S- + S-⊗S+)/2`. Its eigenvalues are
`-2,-1,1`. The pi twist negates only exchange entries. In a half of length `l=n/2`, set

```
B = u h + v h²
A_l = H_open,l + B_12 + B_(l-1,l)
H0 = A_l ⊗ I + I ⊗ A_l + 2 c I
R = h_23 - B_12 - B_34 - c I.
```

If `R >= 0`, the two embedded copies of R at the cuts sum exactly to `H-H0`.
They act on four sites each. For n=8,10,12 their supports are disjoint, although disjointness
is not needed for positivity. For the twisted cut, conjugate R by the diagonal pi rotation
on sites 3 and 4: the crossing bond twists while B_34 is invariant. The script checks this
integer identity. It does not assume that the twisted ring is SU(2)-invariant.

Every two-site SU(2)-invariant self-adjoint operator is a linear combination of `I,h,h²`;
the scalar part is absorbed into c. Thus this is the full **symmetric two-site boundary
family**, but not the set of all conceivable local comparisons. Four-start direct
minimization optimizes the ground bound. A separate 438-point grid/refinement per length
minimizes `(Sz=0 count, total count)` at the original tail cutoff. This finite search does
not certify global count optimality. It does not explore independent asymmetric boundaries.

The selected count parameters are particularly simple:

| n | u | v | c (exact) |
|---:|---:|---:|---:|
| 8 | -3/10 | -1/10 | -2.0307364 |
| 10 | -1/5 | 0 | -2.2154884 |
| 12 | -1/5 | 1/20 | -2.4266500 |

The four-site matrix decomposes by Sz into blocks of dimensions `19,16,10,4,1` for
nonnegative Sz. Negative sectors are related by spin flip. `certify.jl` and `check.jl`
perform exact Rational{BigInt} LDL checks on these blocks, including singular-pivot
handling, and check the twisted conjugation. The JSON stores all pivots as exact strings.
The ground-optimized parameters also pass these checks.

Optimizing the ground bound is the wrong objective for an anchor count. At n=12 it
improves the ground bound from -18.740550 (bare cuts) to -17.527216, but gives 914 states
below the original cutoff across all Sz, versus 773 for bare cuts. The count-optimized
comparison has a weaker ground bound, -17.690178, and only 475 states below that cutoff.
The corresponding Sz=0 counts are 212, 191, and 121 respectively (see generated tables
for the measured values; the total counts include negative Sz).

### What rho and k mean

Let `e_(m,i)` be the ordered spectrum of H0 in sector m, and let
`K_m(rho)=#{i:e_(m,i)<rho}`. Min-max gives

```
lambda_(K_m+1)(H_m) >= lambda_(K_m+1)(H0_m) >= rho.
```

This is an **upper bound on the number of H levels below rho**. It is not an assertion
that H actually has K_m such levels. The physical count is often much smaller. Both
counts, the anchor available at the physical count, the stronger anchor available at
the comparison count, and the first five ordered energy losses are recorded for every
`(n,g,Sz)` in [results/anchor/sectors.csv](results/anchor/sectors.csv). There are 66 nonnegative-sector
rows; negative Sz gives the identical row. `results/anchor/analysis.json` (regenerable, gitignored) also contains losses
for bare-cut, ground-optimized, and count-optimized comparisons separately.

For the absolute shifted-trace budget epsilon=3e-5, the uniform-dimension cutoff is

```
rho_uniform = -n a + log(3^n/epsilon)/beta,
a=700741/500000, beta=49/4.
```

It equals approximately -9.644247, -12.267846, -14.891444 for n=8,10,12. These are not
exactly E0+2, and they use the *entire* error budget for the tail. In particular, n=8
untwisted has only two Sz=0 levels below the absolute-error cutoff, but five below E0+2.
The exact anchor certificate uses rational cutoffs rounded upward to six decimals and
proves the comparison counts agree on both sides of the half-spectrum enclosures.
The transcendental derivation of this baseline cutoff is only numerical here.

The failure in the proposed one-shot counting inference is elementary. Take
`H=diag(-2,-1,10)`, let V consist of the second and third coordinate vectors, and use
`rho=0 <= lambda_3(H)`. The two Ritz vectors have zero residual, yet the lowest Ritz
value -1 is not a lower bound on the true ground energy -2. A coarse count of two below
rho does not identify the missing eigenvalue. In this study, the upper Ritz levels in
a space of dimension K_m generally lie above the comparison anchor, so the simple
separation hypothesis needed to turn that count into ordered lower bounds fails.

This does not prove that every possible one-shot method fails. A certified lower bound
on the orthogonal complement, or a separate inertia/Birman-Schwinger calculation, could
replace continuation. Such a bound is additional work and does not follow from the
scalar counting anchor. For the count-only route actually proposed, keep homotopy in
the plan as a required component rather than an unspecified fallback.

### A cheaper tail using the comparison spectrum

Suppose continuation has established that all but r ring levels lie above a common
target rho. For i>r, min-max and this count give

```
lambda_i(H) >= max(rho,e_i(H0)).
sum_(i>r) exp[-beta(lambda_i(H)+na)]
  <= sum_i exp[-beta(max(rho,e_i(H0))+na)].
```

Bound the last sum with 32 bins of width 1/8 plus one final bin. Bin zero includes all
comparison energies below `rho+1/8`, charged at rho; later bins are charged at their
lower edges. The final bin charges all remaining states at `rho+4`. Multiplicities
are exact sector/tensor-product counts, not the full dimension at every energy.

The following targets reserve 1.5e-5 for the tail and leave 1.5e-5 for low-spectrum errors:

| n | rho | initial K by Sz=0,1,2,3,4 | total including negative Sz |
|---:|---:|---|---:|
| 8 | -9.925 | 28,22,11,2,0 | 98 |
| 10 | -12.687 | 38,32,15,4,0 | 140 |
| 12 | -15.453 | 49,40,20,6,1 | 183 |

These parameters were selected after the uniform-cutoff optimization; no global joint
optimization of boundary parameters and binned tail is claimed. The exact bin populations
use lower half-chain eigenvalue enclosures, so rounding can only weaken the tail bound.
Each scalar bound uses `exp(-x) <= 1/sum_(j=0)^256 x^j/j!` for positive rational x, rounded
up to denominator 10^24. Julia checks all rational arithmetic and the final budget.
The associated analytic Taylor inequality still needs a Lean lemma. No floating exponential
is accepted by the finite checker. The exact tail totals are just below 1.5e-5.

### Homotopy that was actually tested

Use `H(t)=H0+t(H-H0)`; local PSD gives monotonicity. At a stage retain K vectors and a
strict upper counting barrier `rho_bar < lambda_(K+1)(H(t))`. For trial matrix V form

```
M = V' (H(t)-rho_bar) V
N = ((H(t)-rho_bar)V)' ((H(t)-rho_bar)V).
```

Require M negative definite. This proves there are at least K levels below the barrier;
the propagated anchor proves there are at most K. Apply ordinary Ritz to
`(H-rho_bar)^(-1)` in the image space `(H-rho_bar)V`. If its K generalized Ritz values
`tau_i` are negative, `rho_bar+1/tau_i`, in reversed order, are lower bounds on the first
K eigenvalues. Levels whose lower bounds exceed the target can be removed, giving the
next counting barrier. Monotonicity propagates that barrier to the next t. If the count
reaches zero, monotonicity completes the rest of the path without another eigensolve.

The scripts implement this logic with Float64 generalized eigenproblems and adaptive
steps. Every accepted stage checks negativity; rejected steps are halved. A numerical
Hellmann-Feynman predictor chooses the next step but is not used as a proof. Dropping
levels immediately at the target can collapse the margin and stall at degeneracies;
the implemented rule removes them near the current barrier, preserving a buffer.
The final retained K can therefore exceed the number of levels strictly below the target.
Scalar Lanczos also stalled at some degenerate ring endpoints when asked for exactly K
eigenvectors. Oversampling small trial spaces resolved these stalls; the saved runs include
the stalled steps and subsequent resumes. A clean run with the revised solver can be cheaper.

These are genuine end-to-end **numerical tests of the continuation design**, not certified
continuations. Residual norms and acceptance margins are recorded. A Lean implementation
must replace every floating generalized eigenproblem and definiteness test with checked
integer/rational inequalities, account for Gram conditioning, and prove the counting
indices and strict barriers. In particular, the numerical safety padding is not a proof
of enclosure. The generated tables give stages, peak/final K, and timing per sector;
all K values at every stage are available in the corresponding JSON files.

The uniform-tail Sz=0 tests required 75 stages for n=8 and 73 for n=10 (untwisted).
The length-12 uniform pilot was stopped after two accepted stages: its initial
121-vector eigensolve alone took about 63 seconds. It is retained as an explicitly
incomplete pilot, not included among successful continuation tests.

### Other candidate routes

Three or more blocks were tested with bare cuts and the selected two-site boundary
parameters. The tested partitions include 4+2+2 and 3+3+2 at n=8; 4+3+3 and 4+4+2 at
n=10; and 4+4+4, 6+3+3, 6+4+2, and 3+3+3+3 at n=12. Additional cuts increased the counts.
The best sampled results are in the generated table and all samples are in
`results/anchor/analysis.json` (regenerable, gitignored). This is not a global optimization over multiblock parameters.

Anderson overlapping clusters can give a better scalar ground bound. For example the
five-site open-chain ground energy is approximately -5.830213, giving the scalar bound
`H_n >= (n/4) E_open,5 I`. Its largest Sz block is only 51. However, a scalar comparison
has every one of the 3^n levels below the thermal cutoff, so it supplies no useful low-rank
excited-level count. Keeping overlapping cluster operators instead of replacing them by
scalars produces noncommuting terms: their spectrum is not a tensor convolution. Weighted
overlaps do not remove that missing spectral calculation.

An edge-spin-1/2 compression similarly gives variational upper bounds by itself. A lower
comparison requires a certified discarded-space lower bound and coupling control. No
unproved effective edge model was admitted as an anchor. These observations reject these
routes in their cheap, directly computable forms; they do not exclude stronger constructions.

### Half-chain certification and the remaining Lean work

The maximum half-chain Sz dimensions are 19, 51, and 141 for lengths 4, 5, and 6. A
floating eigensolve proposes Q and diagonal D, then both are rounded to denominator 2^32.
The checker evaluates exactly

```
delta = ||Q Q' - I||_infinity
e = ||A - Q D Q'||_infinity
radius = e + 3 ||D|| delta.
```

Here Q denotes the scaled matrix. When delta<1/2, polar decomposition gives an orthogonal
U with `||Q-U|| <= delta` and `||Q||<2`. Thus
`||A-U D U'|| <= radius`, and Weyl encloses the entire ordered spectrum. Integer products
use Int128 only after a BigInt bound proves that every product and partial sum fits;
norm sums and final comparisons use BigInt rationals. This avoids trusting the floating
eigensolver for eigenvalue completeness or accuracy.

For the selected length-six half, the largest raw enclosure radius is about 2.78e-7.
Stored endpoints are rounded outward to a 1e-6 grid. The warm arithmetic check of the
141-dimensional block took roughly 0.04 seconds in an uncontended Julia run; timings in
later concurrent runs vary. All half-spectrum/tensor counts pass exact checks. Recursive
half-chain certification is therefore not needed for this experiment.

The current JSON stores enclosures and local pivots, but regenerates the integer Q proposals
using Julia/LAPACK. A standalone Lean checker must receive or regenerate these Q matrices
by a trusted deterministic recipe and then check them. For the three selected half-chain
families there are 53,346 Q entries: about 417 KiB if stored as Int64, before intervals and
metadata. Including the ground-optimized families doubles that. The current approximately
200 KiB JSON is not evidence that a few-KiB standalone Lean certificate has been achieved.

No Lean runtime has been measured here. Two dense products for all selected half sectors
require about 11.9 million scalar multiply-add terms. At the supplied 3e4 scalar-loop
iterations/second, a naive port already spends roughly seven minutes on this small part;
packed arithmetic remains necessary. For the continuation, the generated work table gives
`sum 3 K(K+1)/2`, the packed-dot count for a straightforward three-Gram-matrix implementation.
At n=12 a 64-bit-lane vector is about 4.72 million bits. Applying a 2.5-second *budget*
per product to those counts gives hours to days, not 15–60 minutes. This is a cost scenario,
not a lower bound or a measured Lean timing; the supplied “under 2.5 seconds” does not
establish that every product costs 2.5 seconds. Batched Gram products, reuse across stages,
or a stronger one-shot complement certificate could substantially change this budget.

The earlier 4.5–7.5k Lean LoC / 1.5–3 week estimate is not supported for the count-only
route. A planning allowance for this anchor/continuation component is 3.5–6k Lean LoC:
roughly 0.6–1k for sector/local-PSD lifting, 0.8–1.4k for exact spectral enclosures,
0.4–0.7k for tensor counts and bins, 1.2–2k for inverse-Ritz/counting, and 0.5–0.9k for
continuation bookkeeping. Depending on overlap with existing infrastructure, budget
approximately 8–13k total Lean LoC and 4–8 weeks for the full seed route. These are
engineering estimates, not code produced or a guarantee. Replace the 15–60 minute kernel
promise with an explicit benchmark milestone; even an optimistic target of 1–8 hours
requires substantial batching beyond the naive continuation implementation.

### Validation and scope

The six numerical shifted traces agree with the reproduction's tabulated centers to less
than 10^-6, providing a useful check of the bond, normalization, magnetization multiplicities,
and pi twist. Small-sector dense diagonalization, Hamiltonian identities, negative-Sz
symmetry, exact-checker rejection tests, and the missing-level counterexample are covered
by `test/anchor.jl`. Small residuals and trace agreement do not certify numerical counts.

The pinned paper source directory was absent in this checkout. The conventions and trace
centers were read from `models/spin1/src/thermal.jl`, and the pinned upstream identity from
`models/spin1/NOTICE`. No paper theorem has been silently treated as formalized. Nothing
under `lean/` or `models/spin1/results/` was changed.

Open risks are the Lean inverse-Ritz/counting development, strict barrier checks near
degeneracies, conditioning and generation of all ring trial spaces, kernel cost of repeated
Gram/residual calculations, and the lack of a global optimum certificate for boundary
parameters. The exact local positivity and small-half spectra are no longer the main
uncertainties. The useful plan change is to retain the optimized comparison, adopt the
binned tail, and make counting/continuation and its kernel benchmark explicit milestones.

<!-- BEGIN ANCHOR TABLES -->

### Numerical and exact-check results

Ring energies and losses in the next table are NON-CERTIFYING. K counts are for the uniform-cutoff comparison.

| n | g | E0(H) | ground-opt loss | count-opt loss | physical k (m=0) | bare K0 | ground-opt K0 | count-opt K0 |
|---|---|---|---|---|---|---|---|
| 8 | I | -11.336956 | 0.663044 | 1.227759 | 2 | 44 | 63 | 34 |
| 8 | P | -11.156192 | 0.843808 | 1.408523 | 3 | 44 | 63 | 34 |
| 10 | I | -14.094130 | 0.515383 | 0.743774 | 5 | 86 | 115 | 62 |
| 10 | P | -13.982142 | 0.627371 | 0.855762 | 6 | 86 | 115 | 62 |
| 12 | I | -16.869556 | 0.657659 | 0.820622 | 8 | 191 | 212 | 121 |
| 12 | P | -16.797640 | 0.729575 | 0.892538 | 8 | 191 | 212 | 121 |

Exact comparison counts and next-level lower anchors (both twists). Decimals here display the rational values in certificate.json. All unlisted sectors have K=0 at both cutoffs.

| n | Sz | uniform K | uniform next rho | binned K | binned next rho |
|---|---|---|---|---|---|
| 8 | 0 | 34 | -9.643366 | 28 | -9.799542 |
| 8 | 1 | 26 | -9.643366 | 22 | -9.799542 |
| 8 | 2 | 13 | -9.643366 | 11 | -9.799542 |
| 8 | 3 | 2 | -9.643365 | 2 | -9.643365 |
| 10 | 0 | 62 | -12.243159 | 38 | -12.657827 |
| 10 | 1 | 50 | -12.243159 | 32 | -12.657827 |
| 10 | 2 | 27 | -12.243159 | 15 | -12.657827 |
| 10 | 3 | 8 | -12.243159 | 4 | -12.473613 |
| 10 | 4 | 1 | -11.954917 | 0 | -12.473613 |
| 12 | 0 | 121 | -14.890385 | 49 | -15.440805 |
| 12 | 1 | 100 | -14.890385 | 40 | -15.440805 |
| 12 | 2 | 56 | -14.890385 | 20 | -15.373806 |
| 12 | 3 | 18 | -14.890385 | 6 | -15.373806 |
| 12 | 4 | 3 | -14.855058 | 1 | -14.961896 |

Exact binned-tail arithmetic totals; the analytic Taylor and spectral inequalities remain to be formalized in Lean.

| n | target rho | tail upper | budget |
|---|---|---|---|
| 8 | -9.925 | 1.48995173385e-05 | 0.000015 |
| 10 | -12.687 | 1.49126012899e-05 | 0.000015 |
| 12 | -15.453 | 1.49383528391e-05 | 0.000015 |

First five ordered losses in Sz=0 for the count-optimized comparison (NON-CERTIFYING). Full per-sector losses and both definitions of k are in [results/anchor/sectors.csv](results/anchor/sectors.csv).

| n | g | loss 1 | loss 2 | loss 3 | loss 4 | loss 5 |
|---|---|---|---|---|---|---|
| 8 | I | 1.227759 | 1.333551 | 2.480392 | 2.245592 | 2.245592 |
| 8 | P | 1.408523 | 1.306997 | 2.215279 | 2.235966 | 2.235966 |
| 10 | I | 0.743774 | 1.268582 | 2.240853 | 1.952080 | 1.952080 |
| 10 | P | 0.855762 | 1.252706 | 2.053984 | 1.937183 | 1.937183 |
| 12 | I | 0.820622 | 1.114069 | 1.970001 | 1.822924 | 1.822924 |
| 12 | P | 0.892538 | 1.104232 | 1.831977 | 1.811095 | 1.811095 |

Best sampled multiblock comparison at the uniform cutoff (NON-CERTIFYING selection).

| n | partition | parameters | Sz=0 K | total K |
|---|---|---|---|---|
| 8 | 3+3+2 | boundary | 104 | 434 |
| 10 | 4+3+3 | boundary | 207 | 901 |
| 12 | 4+4+4 | boundary | 407 | 1819 |

Binned-tail continuation, all initially nonzero sectors (**NON-CERTIFYING**). Final K is the retained count, which may exceed the actual number below the target. Timings include numerical eigensolves, compilation, and contention; they are not Lean timings. Negative Sz shares the same spectrum. Per-stage k and rho are written to the regenerable (gitignored) results/anchor/homotopy_n*_binned.json.

| n | g | Sz | stages | initial → final K | physical k | Julia seconds | min margin |
|---|---|---|---|---|---|---|---|
| 8 | I | 0 | 28 | 28 → 2 | 2 | 21.0 | 5.7e-05 |
| 8 | I | 1 | 24 | 22 → 1 | 1 | 30.0 | 5.3e-05 |
| 8 | I | 2 | 14 | 11 → 0 | 0 | 0.9 | 0.00015 |
| 8 | I | 3 | 3 | 2 → 0 | 0 | 0.0 | 0.012 |
| 8 | P | 0 | 25 | 28 → 3 | 2 | 48.4 | 0.00027 |
| 8 | P | 1 | 17 | 22 → 2 | 2 | 3.4 | 0.00061 |
| 8 | P | 2 | 13 | 11 → 0 | 0 | 1.7 | 6.2e-05 |
| 8 | P | 3 | 2 | 2 → 0 | 0 | 0.0 | 0.026 |
| 10 | I | 0 | 51 | 38 → 2 | 2 | 54.8 | 4.2e-05 |
| 10 | I | 1 | 48 | 32 → 1 | 1 | 23.3 | 4.2e-05 |
| 10 | I | 2 | 23 | 15 → 0 | 0 | 6.2 | 5.2e-05 |
| 10 | I | 3 | 5 | 4 → 0 | 0 | 0.3 | 0.00088 |
| 10 | P | 0 | 42 | 38 → 3 | 3 | 68.1 | 3.9e-05 |
| 10 | P | 1 | 29 | 32 → 2 | 2 | 41.2 | 5.3e-05 |
| 10 | P | 2 | 21 | 15 → 0 | 0 | 7.3 | 5.6e-05 |
| 10 | P | 3 | 4 | 4 → 0 | 0 | 0.6 | 0.0029 |
| 12 | I | 0 | 73 | 49 → 5 | 5 | 666.2 | 1.2e-05 |
| 12 | I | 1 | 63 | 40 → 4 | 4 | 361.3 | 1.6e-05 |
| 12 | I | 2 | 25 | 20 → 1 | 1 | 56.6 | 0.00015 |
| 12 | I | 3 | 9 | 6 → 0 | 0 | 16.8 | 0.00043 |
| 12 | I | 4 | 2 | 1 → 0 | 0 | 1.6 | 0.016 |
| 12 | P | 0 | 72 | 49 → 5 | 5 | 675.5 | 1.1e-05 |
| 12 | P | 1 | 33 | 40 → 4 | 4 | 222.1 | 3.1e-05 |
| 12 | P | 2 | 24 | 20 → 0 | 0 | 46.7 | 0.00016 |
| 12 | P | 3 | 6 | 6 → 0 | 0 | 16.2 | 0.00024 |
| 12 | P | 4 | 2 | 1 → 0 | 0 | 0.5 | 0.012 |

Naive three-Gram packed-dot work for completed binned runs. Only nonnegative Sz is computed; negative Sz is reused. The final column applies a hypothetical 2.5-second budget per dot, not a measured cost.

| n | g | sector-stages | sum K | pair products | budget hours |
|---|---|---|---|---|---|
| 8 | I | 69 | 747 | 19560 | 13.6 |
| 8 | P | 57 | 659 | 17940 | 12.5 |
| 10 | I | 127 | 1858 | 63882 | 44.4 |
| 10 | P | 96 | 1457 | 50649 | 35.2 |
| 12 | I | 172 | 3200 | 140313 | 97.4 |
| 12 | P | 137 | 2616 | 115320 | 80.1 |

<!-- END ANCHOR TABLES -->

## 4. spin1_kernel: Lean seed feasibility review, Mathlib audit and budget

The long companion to [README.md](README.md): the Mathlib audit, the packed arithmetic and its
trust boundary, the packing and anchor experiments, a mathematical review of the proposed
Peierls/Weyl/counting route, benchmarks, and the remaining implementation budget. Lean files
named below now live in `lean/SeedKernel/` (shared lane primitives: `Packed`, `PackingProofs`,
`Obligations`, `PackedTests`), `lean/SeedKernel/Spin1/` (spin-1 operators, `Optimized`,
`Obligations`, `Counting`, `Tests`, `H6/H8/H10`) and `lean/SeedKernel/Bench/` (benchmark-only
files, library `SeedKernelBench`); the prototype's `Packed.lean`, `Obligations.lean` and
`Tests.lean` were split into a shared and a spin-1 part.

The paper excerpts were fetched from the repository-pinned OpenAI commit
`adc7f1241b42e322a6451854ab7e4b4c146bf78a`, files
`preprints/The-periodic-spin-one-Haldane-gap-September-24-2026/build/`
`{initialization,transfer,thermal72,trials}.tex`
(`scripts/fetch_paper.jl` downloads them into the git-ignored `upstream/` and checks
`results/spin1_kernel/paper-source-sha256.txt`). Their inputs confirm 18 thermal tables, including **odd
lengths**, and the π twist about z.

### 1. Mathlib audit

Paths below are relative to `Mathlib/` in the Mathlib checkout pinned by the repository's
`lean/` project (commit `c55e6e78`, the same pin this prototype was written against).
Declaration names are fully qualified. `results/spin1_kernel/mathlib-audit.txt` records matching source
lines; `scripts/audit_mathlib.jl` (formerly `scripts/spin1_kernel/audit_mathlib.sh`) reproduces the search (re-run in the cleanup with identical
output; the negative search, `results/spin1_kernel/mathlib-negative-search.txt`, has zero hits). A wider
manual search for interlacing, minimax and eigenvalue monotonicity also found nothing; the
only related extras are `Matrix.eigenvalue_mem_ball` (Gershgorin, `LinearAlgebra/Matrix/
Gershgorin.lean`) and `LinearMap.IsSymmetric.card_filter_eigenvalues_eq` (multiplicities).
Negative search is evidence of an API gap, not a logical proof that no equivalent statement
can be derived elsewhere.

| Need | Exact declarations already present | Missing work and incremental Lean LoC estimate |
|---|---|---|
| Hermitian diagonalization and eigenvectors | `Matrix.IsHermitian.spectral_theorem`, `Matrix.IsHermitian.eigenvectorBasis`, `Matrix.IsHermitian.mulVec_eigenvectorBasis`, `Matrix.IsHermitian.eigenvectorUnitary` in `Analysis/Matrix/Spectrum.lean` | Thin matrix/operator wrappers, 80–150 |
| Eigenvalue ordering | `Matrix.IsHermitian.eigenvalues₀`, `Matrix.IsHermitian.eigenvalues₀_antitone`; `LinearMap.IsSymmetric.eigenvalues`, `LinearMap.IsSymmetric.eigenvalues_antitone` in `Analysis/InnerProductSpace/Spectrum.lean` | Ascending-index wrapper and reindexing, 50–100. Existing order is **descending**. `Matrix.IsHermitian.eigenvalues` is reindexed onto the original index type and must not be treated as the ascending list. |
| Trace and spectral invariance | `Matrix.IsHermitian.trace_eq_sum_eigenvalues`, `Matrix.IsHermitian.eigenvalues_eq_eigenvalues_iff`, `Matrix.IsHermitian.roots_charpoly_eq_eigenvalues₀` | Unitary/basis transport wrappers, 100–200; repository already has useful wrappers listed below |
| Rayleigh quotient; extremal eigenvalue existence | `ContinuousLinearMap.rayleighQuotient`, `ContinuousLinearMap.rayleighQuotient_le_norm`, `ContinuousLinearMap.norm_eq_iSup_rayleighQuotient`; `LinearMap.IsSymmetric.hasEigenvalue_iInf_of_finiteDimensional`, `LinearMap.IsSymmetric.hasEigenvalue_iSup_of_finiteDimensional` in `Analysis/InnerProductSpace/Rayleigh.lean` | Ground-state Rayleigh upper-bound wrapper, 100–200. These do not give the full indexed Ritz interlacing theorem. |
| Courant–Fischer / indexed min–max | No ready-made declaration found in the spectral/inner-product/matrix API or wider text search | Dimension-intersection argument and both extrema on k-dimensional subspaces, 600–1000 |
| Weyl inequalities | No Hermitian ordered-eigenvalue perturbation/addition theorem found; Weyl-group hits are unrelated | Operator-order monotonicity and norm perturbation, 250–450 after min–max; full two-index Weyl addition is unnecessary for this task |
| Cauchy interlacing / indexed Rayleigh–Ritz | No ready-made declaration found | Orthogonal compression, ascending indices, 200–350 after min–max |
| PSD / Loewner order | `Matrix.le_iff`, `Matrix.nonneg_iff_posSemidef` with scoped `MatrixOrder` in `Analysis/Matrix/Order.lean`; `Matrix.posSemidef_iff_dotProduct_mulVec`, `Matrix.PosSemidef.dotProduct_mulVec_nonneg` in `LinearAlgebra/Matrix/PosDef.lean`; `Matrix.IsHermitian.posSemidef_iff_eigenvalues_nonneg`, `Matrix.PosSemidef.eigenvalues_nonneg` in `Analysis/Matrix/PosDef.lean` | Exact rational LDL/congruence checker still needed, 350–650; positivity of a matrix and efficient decidable certification are distinct |
| Trace monotonicity | `Matrix.PosSemidef.trace_nonneg`, `Matrix.tracePositiveLinearMap` | Trace of an **antitone scalar spectral function** for A >= B follows from ordered eigenvalues, 100–200. Not available merely from trace positivity, and false for arbitrary f. |
| Jensen | `ConvexOn.map_sum_le`, `ConvexOn.map_centerMass_le` in `Analysis/Convex/Jensen.lean`; `convexOn_exp` in `Analysis/Convex/SpecificFunctions/Basic.lean` | Apply scalar Jensen to spectral weights, normalize rows/columns and sum: Peierls 200–350 |
| Nonorthogonal Peierls / Ritz data | Above spectral, PSD, Jensen building blocks | Gram whitening or frame correction; generalized residual and small-matrix enclosure interface, 500–900 |

Relevant existing repository declarations (read-only):

- `CertifiedGaps.trace_exp_smul_eq_sum_of_isHermitian` in `lean/CertifiedGaps/Ring/Model.lean`.
- `CertifiedGaps.eigenvalues_eq_of_conj`, `CertifiedGaps.piTensor_conj`,
  `CertifiedGaps.conj_bondOp_kronecker`, `CertifiedGaps.piTensor_conjTranspose_mul_self`
  in `lean/CertifiedGaps/Ring/Conjugation.lean`.
- `SpinOneHeisenberg.hamiltonian_isHermitian`, `SpinOneHeisenberg.trace_exp_hamiltonian`,
  `SpinOneHeisenberg.partitionFunction_eq_sum` in `lean/CertifiedGaps/Heisenberg/Hamiltonian.lean`.

A focused spectral foundation is approximately **2.8–4.9k Lean lines**, including overlap
between the rows above. Do not add these audit estimates to the module budget below;
they describe the same work from an API perspective.

### 2. Executable prototype and trust boundary

`Packed.lean` implements:

- little-endian base-3 configurations in a full `3^n`-lane Nat, radix `B = 2^w`;
- digit-pattern masks via finite geometric series, avoiding per-configuration mask generation;
- standard-Sz `H v` with two nonnegative limbs representing signed lanes, including a
  pi-twisted seam;
- coordinatewise signed floor using an aligned bias, with masking after the global shift;
- butterfly reversal of vector lanes in O(log d) big operations;
- one big multiplication for each nonnegative dot product (four for a signed pair);
- Gram matrices and product-state power-filter recipes.

`Optimized.lean` hoists masks out of the filter and hoists reversal out of the k-by-k Gram
loop. It also provides reversal masks built by repeated doubling (shifts and ORs),
avoiding the expensive big divisions in the initial reversal. `FastReverse8.lean` checks
this against the original on a non-palindromic regenerated vector. General correctness
is specified by `FastReverseCorrect` and remains open. This removes measured division work.
It does not solve compact-sector addressing.
`PreparedH8.lean` kernel-checks equality of all three optimized
and original H8 moment integers.

The tested trial is

`v = (n I - H_M)^q e_Neel`,

where `H_M` is the Marshall conjugate of the untwisted even-ring Hamiltonian. This is an
exact integer power filter, not a Chebyshev or tuned thermal certificate. Set
`F = nI-H_M`, `G=vᵀv`, `A=vᵀFv`, `C=(Fv)ᵀFv`. Then the intended energy and residual are

`theta = n - A/G`, `residual² = C/G - (A/G)²`.

The integer checker proves `G>0`, `den*A >= (den*n+num)*G`,
`epsDen²*(C*G-A²) <= epsNum²*G²`, and `A²<=C*G`.
That final inequality prevents a Nat subtraction from silently truncating a negative
residual numerator. The raw exact G/A/C values are in `results/spin1_kernel/exact_witnesses.tsv`.
The Julia generator independently computes them using uncompressed sparse rows and BigInts;
its agreement is a regression check, not part of Lean's trust base.

| n | Degree q | Lane bits | Checked intended Ritz upper endpoint | Checked intended residual upper endpoint |
|---|---:|---:|---:|---:|
| 6 | 24 | 262 | -8.495447 | 0.274967 |
| 8 | 32 | 345 | -11.179504 | 0.279904 |
| 10 | 16 | 222 | -13.548710 | 0.632355 |

**These table entries are arithmetic witnesses, not proved physical eigenvalue bounds.**
Only one trial per size is tested; multiple low levels and a full spectral counting
certificate are not delivered. The residuals are far too loose for the thermal tables.

`PackingProofs.lean` proves `radixPack_head`, `radixPack_tail`, `radixPack_extract`,
`add_noCarry`, and `pack_eq_radixPack`. `Obligations.lean` gives precise propositions
`DotCorrect`, `ReverseCorrect`, `FastReverseCorrect`, `MaskCorrect`, `FloorCorrect`,
`HApplyCorrect`, and `FilterCorrect`.
These are `def ... : Prop`, never `axiom` or theorem-by-placeholder.
They are not assumed by the checkers. The bitwise/matrix interpretations remain unproved.

Important explicit bounds:

- Inputs must have no lanes above d, and each limb must have lanes at most M.
- H accumulation is safe if `3*n*M < B`; the positive filter uses `4*n*M < B`.
- Dot products require **every convolution coefficient**, not just the center, below B.
  `d*M² < B` suffices. The trial modules kernel-check
  `3^n*(4*n)^(2*(q+1)) < 2^w` as a conservative arithmetic bound.
- Signed floor requires `s <= biasBits < w`, negative lanes `<=2^biasBits`, and
  positive lanes plus the bias `< B`. The bias must be divisible by `2^s`.
- With d=531441 and lanes of size `2^55`, a dot needs at least **130 bits**, before
  allowing H-vector growth. The 128-bit dense microbenchmark uses tiny lanes and is safe;
  it is not evidence that 128 bits suffices for production Q=2^55.

`Tests.lean` checks signed H against a separate row formula for n=2,3,4 and both twists;
Marshall filtering for n=2,4; negative rounding `[3,-3,17,-17]/2 -> [1,-2,8,-9]`;
reversal; and a Gram matrix. The n=2 case checks the double bond. A deliberately undersized
radix gives a packed dot of 0 instead of the true 18, exhibiting why carry bounds matter.
`Gram6.lean` additionally checks a 4-by-4 Gram matrix on 729-lane vectors.

`Counting.lean` proves the one-Ritz-value **conditional** counting step in scaled integer
coordinates, with Weyl, separation and anchor as explicit premises. It checks a concrete
counterexample to omitting separation. It does not instantiate a physical spectral theorem.

All proof computations use `decide +kernel`. Printed axiom sets are subsets of
`propext`, `Classical.choice`, `Quot.sound`; no `sorryAx` or native-decide axiom is used.

### 3. Packing and anchor experiments

#### Compact-sector addressing

Counts below are exact Julia enumeration of the Sz=0 sector, not kernel-proved counts.
The packing order is increasing base-3 index restricted to that sector. Each entry counts
**distinct** target-rank minus source-rank differences.

| n | Full dimension | Sz=0 dimension | Distinct compact shifts | Closing-bond shifts | Full vector MiB at 128 bits |
|---|---:|---:|---:|---:|---:|
| 6 | 729 | 141 | 32 | 18 | 0.0111 |
| 8 | 6561 | 1107 | 120 | 94 | 0.1001 |
| 10 | 59049 | 8953 | 614 | 574 | 0.9010 |
| 12 | 531441 | 73789 | 3804 | 3744 | 8.1091 |

At n=12 the per-bond counts are `2,4,6,8,10,12,10,8,6,4,2,3744`.
This refutes the straightforward compact-sector implementation of “about 50 masked shifts.”
It is not a lower bound on every possible sector permutation circuit: a fast cyclic rotation,
bit-permutation network, or a different representation could recover a cheap seam.
Such a construction has to be supplied and benchmarked. Full packing costs 7.202 times as
many lanes at n=12; it was chosen here because its transitions really are fixed shifts.

#### Untuned cut anchor (all floating-point results NON-CERTIFYING)

The always-valid elementary comparison is

`H >= H_open(n/2) tensor I + I tensor H_open(n/2) - 4I`,

because each of the two removed exchange bonds has minimum eigenvalue -2. Removing bonds
without this compensation is not a lower comparison. The eigenvalues of a bond are
-2,-1,1, with multiplicities 1,3,5. Twisting a bond by a unitary preserves this fact.

For `a=700741/500000`, `beta=49/4`, `tau=3e-5`, use the conservative full-dimension cutoff
`rho = -n*a + log(3^n/tau)/beta`. This spends the entire tolerance on the tail;
a real error budget must reserve some for retained levels and scalar enclosures.

| n | H0 ground energy | rho | # H0 levels below rho | # below -na |
|---|---:|---:|---:|---:|
| 6 | -10.000000000000 | -7.020648170750 | 64 | 15 |
| 8 | -13.291502622129 | -9.644246899131 | 162 | 32 |
| 10 | -15.660425045542 | -12.267845627512 | 322 | 64 |
| 12 | -18.740549938849 | -14.891444355892 | 773 | 86 |

A large overcount is not merely a request for more Ritz vectors. The counted-anchor lemma
below requires that the kth Ritz value be below the anchor. We diagonalized all Sz blocks
of H6/H8 in Float64, sorted their full spectra with multiplicities, and compared every k:

| n | max over 1<=k<d of lambda_(k+1)(H0) - lambda_k(H) | Best k | Positive margins |
|---|---:|---:|---:|
| 6 | -1.382576818186 | 1 | 0 |
| 8 | -1.445376754516 | 1 | 0 |

Thus **even ideal Ritz vectors cannot satisfy the separation condition with this H0** in
these two examples. This is a robust numerical diagnostic with margins over one energy unit,
not a certified impossibility theorem. Tuned boundary terms were not optimized or excluded.
A proposed tuning should first demonstrate a positive separated margin, with its exact local
PSD certificate, before a recursive Lean implementation is budgeted.

### 4. Mathematical review of proposed steps 1–4

#### Step 1: Peierls with integer trial families

For orthonormal vectors q_i, scalar spectral Jensen gives

`exp(-beta <q_i,H q_i>) <= <q_i, exp(-beta H) q_i>`.

Summing and using the orthogonal projection <= I proves the proposed lower bound. For a
nonorthogonal family the uncorrected sum is false: H=0 on C² and three copies of one unit
vector give 3 on the proposed lower side but trace 2.

Two valid repairs:

1. If V has full column rank, G=V*V>0 and Q=V G^(-1/2) is orthonormal. Diagonalize
   T=G^(-1/2)(V*HV)G^(-1/2); then `Z >= sum exp(-beta theta_i)` for its Ritz values.
   An integer family costs no Gram penalty when handled through this generalized problem.
2. Normalize columns to w_i and prove `WW* <= c I` (equivalently the nonzero spectrum of
   W*W is <=c). Then `Z >= (1/c) sum exp(-beta <w_i,H w_i>)`. The Gram here is the
   **normalized** Gram, not raw V*V. The frame penalty c can destroy the target accuracy.

A near-identity Gram estimate alone does not justify dividing signed Rayleigh energies by
one common endpoint. Use generalized PSD comparisons or a proven shift to positivity.
Rank and normalization checks are mandatory, even though the generator needs no proof.

#### Step 2: the Weyl/counting lemma, and what the anchor must prove

Write the ordered physical eigenvalues **ascending** as lambda_1,...,lambda_d.
For orthonormal Q with k columns, let T=Q*HQ have Ritz eigenvalues theta_1<=...<=theta_k,
R=HQ-QT, and `||R|| <= eps`. Decompose into trial and complement, and delete the
cross block to get `B = T direct_sum C`. Then `||H-B||=||R||`.

A sufficient certificate is:

- `rho <= lambda_(k+1)(H)` with multiplicity, for k<d;
- **strict separation** `theta_k + eps < rho`;
- correctly ordered enclosures for theta_i and a spectral/operator norm residual bound.

Weyl implies `lambda_(k+1)(B) >= rho-eps`. Since all k Ritz eigenvalues are below
`rho-eps`, C cannot contribute an eigenvalue below `rho-eps`. Therefore the lowest k
values of B are exactly theta_i, and Weyl gives

`theta_i-eps <= lambda_i(H) <= theta_i`, for 1<=i<=k,

where the upper bound is Rayleigh–Ritz. Repeated eigenvalues and threshold equality must be
handled explicitly; strict separation is a convenient safe formulation.

Without separation, `H=diag(-100,0)` and Q=e₂ have theta=0, residual=0 and the valid anchor
`rho=0<=lambda_2(H)`, but lambda_1=-100. Without counting, zero residual can similarly
select any excited invariant subspace. Bounding each column residual separately is not a
bound of eps on the whole block: a Frobenius bound, or an operator-norm Gram bound, is needed.

For integer V, put `G=V*V`, `A=V*HV`, `C=(HV)^*(HV)`. The correct generalized residual test is

`C - A G^(-1) A <= eps² G`, with `G>0`.

Whitening turns this into `R*R<=eps²I`. If approximate Ritz coordinates or approximate
normalization are used, include those errors too. A small-matrix LDL/inertia certificate
can avoid computing square roots inside Lean. None of this supplies the missing anchor.
Lehmann–Goerisch is also not anchor-free; it needs a separately justified spectral count or
complement estimate. Calling homotopy a fallback does not supply a crossing-count proof.

#### Tail and trace monotonicity

With lower enclosures l_i for the retained eigenvalues and
`rho<=lambda_(k+1)(H)`, beta>0,

`Z_shifted <= sum_(i<=k) exp(-beta*(l_i+n*a)) + (d-k)*exp(-beta*(rho+n*a))`.

The shift n*a is essential. Use sector dimensions and multiplicities consistently; a
sectorwise k and rho cannot be inserted into a full-space tail without summing all sectors.
The trace monotonicity requested in the task holds for antitone f, including
`f(t)=exp(-beta*t)`. It is false for unrestricted f (take f(t)=t). No operator-monotonicity
of the exponential is needed: ordered eigenvalue monotonicity suffices after taking trace.

At n=12 the crude tail threshold is -14.891444355892, or shifted energy 1.926339644108.
A ground-energy upper estimate cannot by itself replace an absolute certified rho.
Avoid counting “states near a numerical E0” unless all endpoints and multiplicities are
certified. Scalar exp values may be supplied by rational Taylor/interval certificates, but
the analytic exp enclosure theorem still has to be proved; exp does not disappear from the
mathematical trace bound.

Precision is more severe than the energy window suggests. The n=4 untwisted shifted trace
is approximately 124.885379. A uniform eigenvalue error delta changes its trace bound by
roughly beta*Z*delta, giving delta about `1.96e-8` for a 3e-5 budget, before allocating
other errors. For n=12, Z approximately 1.900547 gives about `1.29e-6`. The observed
0.27–0.63 residuals do not test the needed precision. Filter degree, k, and Gram conditioning
must be chosen from a complete error budget, not just a cutoff near E0+2.

#### Local comparison and recursion

A <=81-by-81 local PSD check is useful only after an exact identity decomposes H-H0 into
embeddings of those PSD operators. Boundary terms assigned to the two halves must cancel
in that identity; overlapping supports and scalar allocations cannot be counted twice.
A half-chain eigenvalue list must include certified multiplicities, all relevant Sz sectors,
and a certified remainder count before pair-sum ordering gives an anchor. Recursive anchors
can accumulate losses and rank inflation. The measured untuned comparison is already too
weak for the required separation at n=6 and 8.

A cheaper development order is to search for a usable tuned comparison numerically, emit
small rational local PSD certificates, and test the separation inequality first. If no such
comparison exists, investigate a verified inertia/Schur-complement scheme using the low
subspace of H0, a complement coercivity bound, or a genuinely counted homotopy. These are
new routes requiring their own cost measurements, not implemented fallbacks here.

#### Step 3: untrusted regeneration

It is correct that V may be generated by any terminating recipe, including arbitrary
rounding: Rayleigh–Ritz does not require proving the recipe approximates eigenvectors.
The checker must still establish G>0, valid arithmetic/carry bounds for its interpretation
of V/HV, its Gram/residual inequalities, and the anchor. Floor errors need not be bounded
for V generation, but cannot be ignored if HV or the reported Gram entries are rounded.
This distinction can simplify the generator substantially.

#### Step 4 and connection to the existing Hamiltonian

The existing `SpinOneHeisenberg.hamiltonian` is Cartesian, complex and **untwisted**. A
base-3 matrix is not definitionally this matrix. Required bridges:

1. Let m_j=s_j-1 and identify index `sum s_j*3^j` with `Fin n -> Fin 3`.
2. Check a local Cartesian-to-weight unitary. Compatible columns for m=-1,0,+1 are
   `(1,-i,0)/sqrt(2)`, `(0,0,1)`, `(-1,-i,0)/sqrt(2)`. With the repository's
   `S_alpha=-i epsilon_(alpha s s')`, these give S^z=diag(-1,0,1) and positive sqrt(2)
   ladder entries. Their product cancels the 1/2 in the exchange to give integer hops +1.
3. Lift this identity sitewise using the existing tensor conjugation lemmas. Restrict to
   n>=2; preserve the n=2 double bond. A blanket n=1 formula is unsafe.
4. Define a twisted Hamiltonian: conjugate one endpoint of the closing exchange by
   `diag(-1,1,-1)`. Only its transverse hops change sign. Prove trace/spectrum invariance
   under the global rotation identifying the paper's P with this z-axis pi twist.
5. Prove Sz conservation, direct-sum decomposition, and multiplicity-preserving spectral
   and trace recombination. Sz=0 alone is not enough for a full thermal trace. Avoiding
   SU(2) theory is possible, but then certify every Sz sector.
6. The positive Marshall filter is valid only for **even untwisted** rings. Source 1
   includes n=5,7,9,11 and pi twists; signed generation or a frustrated seam is needed.
   Source 2's third twist has complex phases; an integer pair representation in 1,omega
   is possible, but its conjugation, Gram and norm form `a²-ab+b²` require extra proofs.

The hApply tests exercise twisted signs in the weight representation. They do not prove
any of these global bridges. The paper's symmetry-resolved spatial-transfer I/O/N moments
are also additional structure beyond the already-proved two-purity transfer inequality.

### 5. Benchmark table and extrapolation

Two serial runs, one fresh Lean process per check: the original PR run (Linux container, Intel
Xeon Platinum 8573C, 16 GiB, `results/spin1_kernel/benchmarks_original_linux.tsv`) and the cleanup re-run
(Apple M5, 32 GiB, `results/spin1_kernel/benchmarks.tsv`, `scripts/spin1_kernel/benchmark.jl`). Seconds / peak RSS GiB.

| Check | file | Linux s | Linux GiB | M5 s | M5 GiB |
|---|---|---:|---:|---:|---:|
| Imports / trivial check | Baseline | 4.205 | 0.798 | 3.24 | 0.60 |
| n=6 dense matvec + dot | Bench6 | 2.265 | 0.800 | 0.52 | 0.61 |
| n=8 dense matvec + dot | Bench8 | 2.500 | 0.800 | 0.52 | 0.64 |
| n=10 dense matvec + dot | Bench10 | 6.457 | 1.118 | 1.33 | 1.03 |
| n=12 dense matvec + dot | Bench12 | 72.996 | 5.548 | 26.38 | 3.70 |
| H6, degree 24, 3 moments | H6 | 3.596 | 0.799 | — | — |
| H8, degree 32, 3 moments | H8 | 6.584 | 1.557 | — | — |
| H10, degree 16, 3 moments | H10 | 21.589 | 4.045 | — | — |
| n=12 matvec + checksum | Matvec12 | 27.960 | 3.761 | 11.34 | 2.51 |
| n=12 reversal + dot | Dot12 | 56.997 | 2.617 | 20.64 | 2.28 |
| n=12 prepared masks, matvec + dot | Prepared12 | 60.279 | 4.694 | 28.34 | 3.41 |
| H8 optimized/original equality | PreparedH8 | 8.111 | 2.352 | 1.43 | 2.29 |
| n=6, four-vector Gram | Gram6 | 2.568 | 0.800 | 1.36 | 0.61 |
| Signed/twisted/floor/Gram cases | Tests | 3.852 | 0.810 | — | — |
| n=12 dot, reversed operand supplied | PlainDot12 | 5.307 | 0.823 | 2.46 | 0.77 |
| n=12 dot with division-free reversal | FastDot12 | 9.957 | 2.856 | 2.96 | 2.57 |
| n=12 prepared matvec + fast dot | Fast12 | 30.369 | 5.105 | 9.63 | 3.56 |
| n=8 fast/original reversal equality | FastReverse8 | 2.785 | 0.800 | 1.24 | 0.67 |
| n=12 prepared matvec + checksum | PreparedMatvec12 | 18.042 | 2.961 | 6.98 | 2.99 |

H6, H8, H10 and Tests are now library modules built by `lake build`; on the M5 a clean build
of the whole library (all of them) takes 8.5 s, with H10 the largest module at 3.2 s.

The dense microbenchmark computes one filter matvec on a full repunit vector, reverses a
vector, multiplies, and checks the closed sum `17*n*3^(n-2)`. Every theorem is elaborated
and kernel-checked afresh. The Hn tests regenerate the entire trial and three dot products.
Times include Lean/Lake startup and imports, elaboration, proof construction and checking;
there is no claim that the table isolates kernel normalization time. RSS is the maximum
reported by Linux child rusage; startup contributes a large floor. Small-n timing variance
is dominated by startup. The final serial run avoids concurrent eigensolvers or bignum jobs.

One 4-million-bit multiplication does not characterize a whole matvec. At n=12 this
prototype uses 68,024,448-bit full vectors at 128 bits per lane, repeated divisions for
masks, many shifts/ands, and reversal. Peak resident memory is orders of magnitude larger
than one packed vector. Hoisting masks and reversals is necessary. The optimized
n=12 dense matvec-plus-dot takes **30.369 s** versus **72.996 s** initially
(2.40x faster), while peak RSS is **5.105 GiB** versus **5.548 GiB**. A dot with a
pre-reversed operand takes **5.307 s**; the division-free reversal-plus-dot takes **9.957 s**,
versus **56.997 s** for the first reversal. Hoisting reversed vectors out of a Gram loop
is consequently a material optimization, not a cosmetic one.

Single-check peak memory
cannot be extrapolated linearly to thousands of theorems; checkpointed modules may be needed.

For illustrative full-packing workloads with k=50 and q=200 (Linux timings; the M5 re-run
is 2.8–3.2x faster at n=12, which divides the hours below but does not change the conclusion):

- `q*k = 10000` matvecs. At roughly 30 seconds per n=12 matvec, this is about **83 hours**
  for one n/twist, before Gram work. This is a conditional workload calculation, not a
  measured production runtime; q and k are not established by these experiments.
- The prepared-mask matvec plus checksum measures **18.042 seconds**, including mask
  preparation. Repeating it 10000 times would take **50.1 hours**. Reusing masks across
  all filter iterations can reduce that cost, but was not measured at n=12 over q=200.
- Three symmetric Gram families need up to `3*k*(k+1)/2 = 3825` dot products. At the
  observed 5.307 seconds per dot with an already-reversed operand, a naive individual-check
  extrapolation is **5.64 hours**, before reversal/preparation. Startup amortization and
  shared expressions can change this; it is not a measured batched Gram runtime.
- Using the optimized full matvec-plus-dot as a 30.369-second budget gives **84.4 hours**
  for 10000 such operations. A complete filter should reuse prepared masks and omit dots
  at intermediate steps, so this is not a lower bound. There is no measured 200-step,
  50-vector n=12 kernel run proving that it fits memory or meets this estimate.
- Compact sector packing, reusable masks and reversals, narrower certified lanes, and
  carefully selected k/q could reduce this substantially. The sector seam circuit and
  anchor are currently unknown, so a defensible optimized total is **unknown**, not a
  confident time interval. Redesigning n_max=10 or 8 is strongly favored by the measured
  size growth, but must revalidate the full analytic seed and precision budget.

### 6. Module-by-module remaining implementation budget

These are engineering estimates for reviewed Lean source, not measured coding throughput.
Ranges account for missing reusable APIs and proof representation choices. They exclude
this small prototype and do not include a general-purpose full two-index Weyl theorem.

| Proposed module | New Lean LoC | Why it is needed |
|---|---:|---|
| PackedRadix / Carry / SignedFloor | 350–650 | Connect proved reference extraction to shift/mask operations, convolution carry bounds |
| Masks / Reverse / Gram | 350–650 | Geometric masks, reversal, signed products, reuse and Gram interpretation |
| SpinOnePacked / SectorDecomposition | 600–1000 | H action, seams, bounded accumulation, all Sz sectors |
| CartesianWeightBridge | 400–700 | Local sqrt(2) unitary and global indexed matrix transport |
| OrderedMinMax | 700–1100 | Ascending eigenvalues, subspace dimension arguments, Ritz |
| EigenvalueOrder / Weyl / Interlacing | 450–750 | Loewner monotonicity, operator perturbation, compression |
| RationalPSD / GeneralizedRitz | 500–850 | LDL/congruence certificates, G>0, rational spectral intervals |
| ResidualCounting | 300–500 | Generalized block residual and separated counted-anchor theorem |
| Peierls / TraceTail / ScalarExp | 500–900 | Nonorthogonal treatment, finite trace/tail and exp enclosures |
| LocalAnchor / RecursiveCount | 900–1800 | Local comparison identity, pair sums, all-sector recursive certification |
| Recipe / CertificateInterface | 250–450 | Bound checks, failure cases, replay and theorem-facing interface |
| ThermalTableAssembly | 300–600 | 18 enclosures, shifts, multiplicities, error budgeting |
| **Thermal route subtotal** | **5600–9950** | Does not solve the numerical anchor-design problem |
| PeriodicMPS / ExactTransfer | 1000–1800 | Norm and insertion contractions, integer matrix/block checker |
| SymmetryTransfer / CappedMass / DominantIsolation | 1600–3000 | I/O/N moment identities and finite-seed analytic bridge |
| SeedAssembly / ExistingChainGlue | 250–450 | Feed (36,49/4) seed into existing doubling chain |
| **Whole seed-route total** | **8450–15200** | Source-1 scope; extra third-twist work for source 2 |

No 1.5–3-week completion promise is justified. Conditional on finding a usable anchor and
representation, a planning range is **8–16 expert weeks** for this broader proof work,
plus potentially unbounded mathematical/certificate search. The estimate corresponds to
roughly 200 reviewed Lean lines per working day across the 8450–15200-line range;
that throughput is an explicit planning assumption, not a measured productivity result. This is
a staffing estimate, not evidence that an anchor exists.
A first milestone should be a fully bridged H6
counted enclosure, before expanding to H8 or implementing all transfer/trial modules.

### 7. Plan changes and open risks

1. Put the anchor/separation test first. The untuned comparator fails before precision
   becomes relevant. Require a concrete tuned local certificate and counted margin.
2. Replace the assumed compact O(50)-shift matvec with an actual seam circuit, or budget
   the measured full-base-3 costs. Cache mask construction and vector reversal.
3. Use generalized Ritz eigenvalues for the lower trace; avoid an uncontrolled frame
   penalty. Budget residual precision against beta times the shifted trace.
4. Prove a reusable ordered min–max/Weyl layer and the Cartesian/twist bridge. Do not
   present arithmetic witnesses as physical eigenvalue theorems.
5. Make the first complete acceptance test H6: all sectors, a counted lower enclosure,
   a Ritz upper enclosure, a shifted trace interval, and `#print axioms` on that theorem.
6. Explore n_max<=10 or source 2 only after charging its complex twist, trial length,
   analytic transfer identities and precision. The existing optimized scalar note is
   useful evidence for redesign, not a Lean-certified shortcut by itself.

Remaining risks: tuned anchor existence, rank inflation, norm-vs-column residual errors,
Gram conditioning, odd rings and twists, degeneracies and threshold equalities, all-sector
multiplicities, carry/borrow proofs, expression-retention memory, Mathlib API integration,
scalar exp bounds, and the unimplemented analytic/MPS seed bridge. The project deliberately
keeps all of these visible; none is hidden behind an extra axiom or an unchecked generator.

### 8. Validation record

Re-run in the cleanup of the standalone prototype (Apple M5, Julia 1.12, Lean v4.35.0-rc3); the
consolidated checks and their times are in README "Reproduce".

- The standalone Std-only Lake project built from clean in 8.5 s (the spin-1 modules of
  `SeedKernel` still import no Mathlib).
- `#print axioms` for all 11 audited declarations (now `lean/SeedKernel/Audit.lean`; also
  `results/spin1_kernel/axioms.txt`) shows only subsets of `propext`, `Classical.choice`, `Quot.sound`;
  no `sorryAx`, no `Lean.ofReduceBool` (no `native_decide`).
- Every `lean/SeedKernel/Bench` file checks (exit 0), the slowest being Bench12 (26 s) and Prepared12
  (28 s); none needed the 5-minute limit.
- The generator (now `src/spin1_kernel.jl`, run by `certify.jl spin1_kernel`) regenerates
  `lean/SeedKernel/Spin1/H{6,8,10}.lean` and
  `results/spin1_kernel/exact_witnesses.tsv` byte-identically (105 s). The witnesses are tight: raising
  the H6 Ritz numerator by 1 (8495447 to 8495448) makes `decide +kernel` report the
  proposition false.
- `scripts/spin1_kernel/experiments.jl` reproduces the structure and anchor tables exactly and the
  Float64 level table to about 1e-13.
- `test/spin1_kernel.jl` (24 exact tests, 2 s) checks sparse-row/packed agreement, digit masks,
  reversal and dots for n = 2..6, the carry-overflow counterexample, and the H6 recipe.
- `check.jl` verifies the source hashes in `results/spin1_kernel/witnesses.json` and the Lean
  witness files against the stored witnesses, builds `SeedKernel` and checks the axiom audit;
  `--full` regenerates the three witnesses, rebuilds `SeedKernel` from clean and checks every
  benchmark file.

The float eigenvalue/anchor tables remain non-certifying.

## 5. ising: the two-level route, the all-model inventory and the remaining proof work

This is the long companion to [README.md](README.md). Everything below that is not a Lean
theorem is labelled; the floating-point numbers are planning evidence only.

### A cheaper spectral route for the selected point

The selected point is the longitudinal Ising ring `H = −Σ ZZ − g Σ X − h Σ Z` at
`g = 1/4, h = 1/20`, seed `(n₀, β₀) = (6, 97/4)`. Two approximate eigenvectors per ring plus a
classical diagonal comparison suffice numerically; full Rayleigh–Ritz matrices, a recursive
half-chain anchor, momentum sectors and SU(2) are unnecessary at this point.

**Classical anchor.** Let `D` be the classical diagonal part of `H`. Since `−gX ≥ −gI` on one
spin, `H ≥ D − ngI`. Exactly two periodic configurations (all up, all down) have no domain
wall; every other one has at least two, so its classical energy is at least
`−n(1+h) + 4 + 2h`. At `g = 1/4, h = 1/20`:

```
ρₙ = −n(1+h+g) + 4 + 2h = (82 − 26n)/20,     λ₃(H₆) ≥ −3.7,     λ₃(H₁₂) ≥ −11.5.
```

The direct classical counts (`anchor6`, `anchor12`) are proved in Lean. The operator
inequality and the eigenvalue-count consequence (min–max) are still proof obligations.

**Residual intervals.** For the integer Hamiltonian `A = 20H` and scale `S = 2⁵²`, the kernel
checks

```
|vᵢᵀvᵢ − S²| ≤ 10⁻¹² S²,   |v₁ᵀv₂| ≤ 10⁻¹² S²,   Σᵢ ‖S A vᵢ − λᵢ vᵢ‖² ≤ 10⁻¹⁸ S⁴,
λ₁ + S < λ₂,               λ₂ + S < (20ρ) S.
```

Once the norm identities are proved, each residual gives an eigenvalue within
`10⁻⁹/(20√(1−10⁻¹²)) < 5.01×10⁻¹¹` of the candidate `μᵢ = λᵢ/(20S)`. The two intervals are
disjoint and lie below ρ, so a proved count of at most two eigenvalues below ρ makes them
exactly the first two. Cross-Gram orthogonality is then unnecessary; it is checked only for
comparison. For larger spectra the analogous saving is to certify small clusters of
overlapping intervals and their multiplicities instead of a full k×k Gram matrix; that still
needs a counting theorem and a reliable anchor.

**Trace sandwich.** With a deliberately conservative energy radius `r = 10⁻¹⁰`,

```
Σᵢ₌₁² exp(−β(μᵢ + r)) ≤ Zₙ(β) ≤ Σᵢ₌₁² exp(−β(μᵢ − r)) + (2ⁿ − 2) exp(−βρₙ).
```

`scripts/ising/feasibility.jl` evaluates this in Float64 (**non-certifying**; it assumes the spectral
theorem above), writing `results/ising/feasibility.json`:

| Quantity | Sandwich result | Required target |
|---|---:|---:|
| Spatial defect upper bound | 1.210665×10⁻⁶ | 9.064×10⁻⁶ |
| Thermal defect upper bound | 9.949076×10⁻⁹ | 3.141×10⁻⁵ |
| 6-site normalised tail at β | 2.7771×10⁻²⁷ | included above |
| 12-site normalised tail at β | 1.2424×10⁻¹⁰ | included above |

The predicted log-partition intervals also fit strictly inside the four endpoints that
`LogEnclosures` asks for, so what remains is a proof task, not a need for sharper vectors.
The `h = 1/10` point has only one level below the inventory cutoff, but the rank-two classical
anchor does not isolate that level by itself; `h = 1/20` is the cheapest point found with a
complete numerical anchor argument.

**Scalar exponentials.** A concrete route is `Real.exp_bound` on each exponent divided by 1024
(the largest exponent here is about 620), followed by ten outward-rounded rational squarings.
A 16-term Taylor proposal is amply sufficient at these margins, but its rounding and remainder
checks must themselves be proved.

### Inventory of all 43 seed targets

`scripts/ising/inventory.jl` reads the stored spectra of every certified model (30 longitudinal
Ising points, 5 staggered-Heisenberg, 3 spin-½ dimerised, 2 large-D, and one each of the
Heisenberg ladder, spin-1 dimerised δ = 9/10 and Ising ladder g = 3) and writes
`results/ising/inventory.{json,tsv}`. The seed defects are not eigenvalue tolerances; to avoid
silently using all of a seed's slack it sets

```
σ = min(εₛ − εₛ,source, εₜ − εₜ,source),   η = min(10⁻⁶, σ/16),   r = η/(4β),
```

and counts the levels below `E₀ + (log D + log(1/η))/β`. "source" is the stored external
enclosure, not a Lean theorem. Counts use the nominal stored spectra with multiplicities and
are **non-certifying**.

| Family | Largest full dimension | Largest Sz block, when used | k at the cutoff | Full three-Gram arithmetic scenario |
|---|---:|---:|---:|---:|
| Longitudinal Ising, 30 points | 4,096 | 4,096 (no Sz reduction) | 1–145 | 0.014–150 s per point |
| Spin-1 dimerised δ = 0.9 | 6,561 | 1,107 | 195 | 5–80 s |
| Staggered Heisenberg, 5 points | 65,536 | 12,870 | 424–1,249 | 0.14–19.8 h |
| Heisenberg ladder J = 4 | 65,536 | 12,870 | 1,238 | 1.16–19.4 h |
| Spin-½ dimerised, 3 points | 65,536 | 12,870 | 862–1,558 | 0.56–30.8 h |
| Ising ladder g = 3 | 65,536 | 65,536 (Sz not conserved) | 322 | 15.7 h |
| Large D = 4 | 531,441 | 73,789 | 1,480 | 15.4–379 h |
| Large D = 3.206 | 531,441 | 73,789 | 2,465 | 42.5–1,051 h |

These are **arithmetic-only scenario estimates** (`scripts/ising/costs.jl`), not measured certificate
run times or proven lower bounds. For conserved charge they range from k spread evenly over
the charge sectors to all k in the largest block, pricing each product at the largest block
width. They exclude vector generation, orthogonalisation, packing, anchors, spectral proofs
and scalar exponential bounds; no momentum, cyclotomic or SU(2) saving is assumed. The low β
of the large-D seeds is why thousands of levels matter despite small seed lengths. These
numbers reject a blanket 15–60 minute estimate for these models with full k×k Gram checks;
they do not refute a separate spin-1 Heisenberg proposal at β = 49/4.

#### Kernel multiplication benchmark

`scripts/ising/costs.jl` checks `(2ᵇ−1)(2ᵇ−3) = 2²ᵇ − 2ᵇ⁺² + 3` by `decide +kernel` on dense
operands (Linux run recorded in `results/ising/costs.json`):

| Operand bits | Kernel time |
|---:|---:|
| 262,144 | 0.00111 s |
| 1,179,648 | 0.00471 s |
| 3,706,560 | 0.0304 s |
| 21,251,232 | 0.415 s |

These are primitive microbenchmarks; `costs.jl` extrapolates with a local size exponent 1.15
from the nearest measured size.

#### Full table

`julia --project=. scripts/ising/report_tables.jl` prints this table from `results/ising/costs.json`.
`η` and `r` as above; Gram minutes exclude vector generation, anchors, analytic proofs and
encoding costs.

| Model / parameter | D | β | εₛ | εₜ | η | r | cutoff − E₀ | k | Gram minutes, balanced–concentrated |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| LargeD 4_1 | 531441 | 2.25 | 0.06411 | 0.04053 | 3.79e-08 | 4.21e-09 | 13.454 | 1480 | 924–2.27e+04 |
| LargeD 1603_500 | 531441 | 3.048 | 0.19164415 | 0.046462453 | 3.66e-11 | 3e-12 | 12.21 | 2465 | 2.55e+03–6.31e+04 |
| Heisenberg ladder J4 | 65536 | 2 | 0.03509 | 0.02463 | 3.63e-07 | 4.53e-08 | 12.96 | 1238 | 69.5–1.17e+03 |
| Staggered h1_1 | 65536 | 5.25 | 0.0006937 | 0.0007327 | 8.26e-10 | 3.94e-11 | 6.096 | 424 | 8.36–137 |
| Staggered h3_4 | 65536 | 5.25 | 0.002343 | 0.003123 | 3e-08 | 1.43e-09 | 5.4117 | 534 | 13.2–217 |
| Staggered h1_2 | 65536 | 5.5 | 0.01075 | 0.01116 | 1.87e-07 | 8.51e-09 | 4.8331 | 786 | 28.2–470 |
| Staggered h7_20 | 65536 | 5.75 | 0.03266 | 0.02744 | 2.56e-07 | 1.11e-08 | 4.5686 | 1152 | 60.2–1.01e+03 |
| Staggered h1_4 | 65536 | 6.25 | 0.07987 | 0.04477 | 2.22e-07 | 8.87e-09 | 4.2259 | 1249 | 70.7–1.19e+03 |
| Dimerized half d1_2 | 65536 | 6 | 0.01133 | 0.01224 | 1.94e-07 | 8.07e-09 | 4.4247 | 862 | 33.9–565 |
| Dimerized half d3_10 | 65536 | 6.5 | 0.09516 | 0.0518 | 5.3e-08 | 2.04e-09 | 4.2837 | 1558 | 110–1.85e+03 |
| Dimerized half d1_4 | 65536 | 7.25 | 0.1658 | 0.05123 | 1.32e-07 | 4.55e-09 | 3.7146 | 913 | 38–634 |
| Dimerized one d9/10 | 6561 | 3.5 | 0.04972 | 0.03206 | 2.71e-07 | 1.94e-08 | 6.8315 | 195 | 0.0845–1.33 |
| Ising ladder G3 | 65536 | 1.5 | 0.1125 | 0.05162 | 2.98e-07 | 4.96e-08 | 17.412 | 322 | 941–941 |
| Ising g1/4 h1/20 | 4096 | 24.25 | 9.064e-06 | 3.141e-05 | 1.41e-11 | 1.45e-13 | 1.3734 | 2 | 0.000706–0.000706 |
| Ising g1/4 h1/10 | 4096 | 13.75 | 4.675e-06 | 1.799e-05 | 3.47e-11 | 6.31e-13 | 2.3565 | 1 | 0.000235–0.000235 |
| Ising g1/4 h1/4 | 4096 | 6.5 | 2.226e-06 | 8.857e-06 | 1.95e-11 | 7.5e-13 | 5.0736 | 25 | 0.0765–0.0765 |
| Ising g1/4 h1/2 | 4096 | 5.5 | 2.063e-06 | 8.242e-06 | 4.56e-12 | 2.07e-13 | 6.2602 | 25 | 0.0765–0.0765 |
| Ising g1/4 h1/1 | 4096 | 5.25 | 2.494e-06 | 9.964e-06 | 3.56e-12 | 1.7e-13 | 6.6053 | 13 | 0.0214–0.0214 |
| Ising g1/2 h1/20 | 4096 | 13 | 0.001314 | 1.872e-05 | 3.82e-12 | 7.35e-14 | 2.6622 | 3 | 0.00141–0.00141 |
| Ising g1/2 h1/10 | 4096 | 8 | 0.0002267 | 1.129e-05 | 5.98e-10 | 1.87e-11 | 3.6945 | 12 | 0.0184–0.0184 |
| Ising g1/2 h1/4 | 4096 | 4.5 | 7.237e-06 | 7.138e-06 | 3.25e-11 | 1.8e-12 | 7.2153 | 74 | 0.654–0.654 |
| Ising g1/2 h1/2 | 4096 | 4 | 1.625e-06 | 6.517e-06 | 1.48e-12 | 9.24e-14 | 8.8894 | 49 | 0.288–0.288 |
| Ising g1/2 h1/1 | 4096 | 4 | 1.974e-06 | 7.885e-06 | 7.46e-12 | 4.66e-13 | 8.4847 | 25 | 0.0765–0.0765 |
| Ising g3/4 h1/20 | 4096 | 8 | 0.0371 | 0.0003654 | 4.67e-09 | 1.46e-10 | 3.4374 | 21 | 0.0544–0.0544 |
| Ising g3/4 h1/10 | 4096 | 5.25 | 0.008573 | 6.931e-05 | 2.96e-10 | 1.41e-11 | 5.7632 | 72 | 0.619–0.619 |
| Ising g3/4 h1/4 | 4096 | 3.5 | 0.0002372 | 0.0001331 | 3.16e-09 | 2.26e-10 | 7.9685 | 103 | 1.26–1.26 |
| Ising g3/4 h1/2 | 4096 | 3.5 | 4.348e-06 | 1.065e-05 | 4.52e-11 | 3.23e-12 | 9.1823 | 68 | 0.552–0.552 |
| Ising g3/4 h1/1 | 4096 | 3.5 | 1.856e-06 | 7.342e-06 | 2.56e-11 | 1.83e-12 | 9.3445 | 25 | 0.0765–0.0765 |
| Ising g1/1 h1/20 | 4096 | 6.5 | 0.1841 | 0.007425 | 5.97e-08 | 2.3e-09 | 3.8387 | 25 | 0.0765–0.0765 |
| Ising g1/1 h1/10 | 4096 | 4.25 | 0.05123 | 0.004667 | 2.98e-08 | 1.75e-09 | 6.0344 | 75 | 0.671–0.671 |
| Ising g1/1 h1/4 | 4096 | 3 | 0.002469 | 0.002245 | 2.01e-08 | 1.68e-09 | 8.68 | 145 | 2.49–2.49 |
| Ising g1/1 h1/2 | 4096 | 3 | 5.079e-05 | 9.929e-05 | 1.1e-10 | 9.17e-12 | 10.416 | 122 | 1.77–1.77 |
| Ising g1/1 h1/1 | 4096 | 3 | 2.207e-06 | 7.69e-06 | 5.42e-11 | 4.51e-12 | 10.652 | 37 | 0.166–0.166 |
| Ising g3/2 h1/20 | 4096 | 3.5 | 0.1478 | 0.05424 | 4.36e-07 | 3.11e-08 | 6.5612 | 56 | 0.376–0.376 |
| Ising g3/2 h1/10 | 4096 | 3.25 | 0.09848 | 0.04085 | 5.67e-07 | 4.36e-08 | 6.985 | 60 | 0.431–0.431 |
| Ising g3/2 h1/4 | 4096 | 2.75 | 0.01818 | 0.01278 | 1.55e-07 | 1.41e-08 | 8.7271 | 92 | 1.01–1.01 |
| Ising g3/2 h1/2 | 4096 | 2.5 | 0.001195 | 0.001673 | 1.52e-08 | 1.52e-09 | 10.529 | 93 | 1.03–1.03 |
| Ising g3/2 h1/1 | 4096 | 2.5 | 2.241e-05 | 3.345e-05 | 4.57e-10 | 4.57e-11 | 11.93 | 61 | 0.445–0.445 |
| Ising g2/1 h1/20 | 4096 | 2.25 | 0.03252 | 0.03673 | 1.73e-07 | 1.92e-08 | 10.618 | 110 | 1.44–1.44 |
| Ising g2/1 h1/10 | 4096 | 2.25 | 0.02941 | 0.03285 | 2.93e-07 | 3.25e-08 | 10.383 | 92 | 1.01–1.01 |
| Ising g2/1 h1/4 | 4096 | 2.25 | 0.01578 | 0.01669 | 3.22e-07 | 3.58e-08 | 10.341 | 81 | 0.782–0.782 |
| Ising g2/1 h1/2 | 4096 | 2.25 | 0.003329 | 0.003196 | 1.57e-08 | 1.75e-09 | 11.683 | 85 | 0.861–0.861 |
| Ising g2/1 h1/1 | 4096 | 2 | 0.0001302 | 0.0003862 | 2.91e-09 | 3.64e-10 | 13.986 | 79 | 0.744–0.744 |

### Implementation observations

- A single-expression direct 12-site checker was stopped at 7.4 GB RSS and 38 CPU s without
  completing. Splitting the scalar work into 64 independently checked blocks
  (`GeneratedChunks`) succeeds; `ChunkCoverage` proves the block list is exactly the fixed
  range enumeration, so a generator cannot duplicate or omit a block.
- The first packed version checked the vector arithmetic in about 0.4 s but spent about 140 s
  elaborating two padded decimal literals. Storing compact 56-bit slots and widening them to
  288-bit lanes inside Lean (`Packed.widen`) cut elaboration to about 6 s.
- Do not treat packed primitive timings as full-checker timings: vector preparation, Gram
  products, classical counts, scalar bounds and peak kernel memory are separate costs.

### Remaining proof work for the selected point

The missing theorem is `LogEnclosures` (README, "Seeds and certified bounds"). A route to it:

| Component | Additional Lean lines (estimate) | Obligation |
|---|---:|---|
| Matrix identification and scalar-sum interpretation | 200–400 | Binary states ↔ `Fin n → Fin 2`; prove `A = 20H` against `RPRing.hamiltonian`; identify the checked sums as norms and Gram entries |
| Residual-to-spectrum bridge | 150–350 | Apply `exists_near_of_weighted_residual` through the Hermitian spectral theorem |
| Classical comparison and counting | 250–500 | Prove `H ≥ D − ngI`, transfer the classical count by min–max, exhaust two eigenvalues |
| Rational exponential and trace bounds | 200–450 | Outward scalar exp bounds and the trace sandwich; derive the four log bounds |
| Integration | 100–250 | Connect to `IsingLongitudinal.gap` without the seed `sorry` |
| Optional packed-backend correctness | 300–700 | Mask/flip/reverse, carries, middle-product extraction, widening; unnecessary with the chunk route |

That is an estimated 900–1,950 lines by the direct route (1,200–2,650 with a reusable packed
backend), not a measured completion. A backend covering all seven model families needs
symmetry identification, anchors, cluster multiplicities and scalar bounds; 8–15k lines is the
more defensible planning figure. The larger Gram scenarios rule out a uniform minute-scale
run time without further structural reductions.

### Suggested order

1. The selected longitudinal Ising point: explicit two-state anchor and a large seed margin.
   A separate theorem for this point can become unconditional before the 30-point `seeds`.
2. The other longitudinal Ising points and spin-1 dimerisation δ = 0.9 (largest Sz block 1,107).
   General points still need anchors; the spin-1 Sz basis change is extra work.
3. Strong staggered fields (h = 1, 3/4), then the Heisenberg ladder and the easier spin-½
   dimerisation targets, prioritised by cluster counts and anchor strength.
4. Weaker staggered fields, harder dimerisation targets and the Ising ladder (no Sz
   decomposition; the most expensive Gram scenario among the spin-½ models).
5. Large D = 4, then D = 3.206; the latter has only about 5.85×10⁻¹⁰ slack against its stored
   enclosure, so redesigning its seed is more promising than a dense Gram plan.
