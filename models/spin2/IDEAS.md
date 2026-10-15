# Spin-2: reassessing the thermal starting point

## Reproducible audit for review, 2026-10-14

The implementation and October 11 strategy memos are already on `main`. The remaining
reviewable contribution is a standalone Julia replay of the two numerical tables below.
Reusing the existing scalar screen is cheapest, but it writes only the polynomial table;
replaying the full MPO trajectory would test a different claim at much higher cost.
Choose a bounded audit of the scalar majorant and the four-site Gibbs/Trotter comparison.
The polynomial scan is scalar in n and beta, with cost controlled by degree, precision,
and the scanned powers. The four-site reference enumerates fixed magnetization blocks;
it is deliberately not a method for large rings. No new MPO construction is claimed.

## Big-picture decision, 2026-10-11

The objective is one admissible, independently checked finite seed for the uniform-gap
criterion. The immediate deliverable is a measured feasibility decision for producing that
seed, not another long high-temperature trajectory.

| Route | Reason to consider it | Scaling risk and decision |
|---|---|---|
| Local gates, final-trace certification, and charge blocks | Reuses implemented rigorous ingredients and avoids the extensive initial Taylor restriction | At fixed retained ranks, work grows with chain length and the number of product steps; ranks and certified residuals may grow with both length and beta. Chosen first, subject to measured error and verification budgets |
| Joint certification of purity ratios or spatial transfer contraction | Could exploit cancellation of the extensive free energy and reduce dependence on ring length | A new correlated error or contraction theorem is required; a compressed transfer spectrum alone does not certify the exact quantum model. Develop only with an explicit checkable remainder |
| More accurate trial states or symmetry alone | Improves the energy shift or cost of existing contractions | Neither controls the unaccounted thermal states; support the selected route rather than serve as the main proof strategy |

First map promising seed sizes and temperatures with exploratory Julia/ITensors calculations
and evaluate the rigorous criterion conditionally. Then test the chosen certificate at
meaningful temperatures on small exact references and increasing MPO lengths. Track total
trace overestimation, rank, memory, and independent checker time separately. Calibrate the
same construction on spin 1 using its existing certificate as a reference; success there is
necessary evidence of practicality, not evidence that spin 2 will follow automatically.

For log partition-function errors, spatial cancellation is potentially useful:
`log(P_s) = log(Z_2n) - 2log(Z_n)`. Independent interval errors still add, however, and
thermal purity compares different temperatures, so its extensive free-energy terms do not
cancel in the same way. Any proposed cancellation must be proved in the certificate.

Proceed to a large seed only after the combined splitting and compression budget and the
observed rank growth support it. If they do not, record the failure and prioritize a
transfer/contraction or correlated-ratio bound instead of extending the quadratic run.
No spin-2 gap certificate is claimed by this strategy memo.

Date: 2026-10-11. Scope: strategy and small rigorous scalar diagnostics; no new large-ring
contractions. Existing `results/` files remain read-only. No spin-2 gap is certified.

The current 256-site quadratic-polynomial route starts at a temperature scale near 10^-7
and requests 28 squarings. Before extending it, separate restrictions imposed by the
polynomial approximation from losses introduced by MPO compression and its error bound.

## Candidate routes

| Route | Why it might help | Obstruction and cost |
|---|---|---|
| Higher-degree initial polynomial with adaptive scaling | Spend more effort once to remove many near-identity squarings; the existing scalar remainder can be checked cheaply | A good scalar approximation does not imply a cheap MPO. Explicit powers of H can grow exponentially in degree; compressed construction needs its own rigorous error |
| Local even/odd bond product with a trace upper-bound inequality | Starting gates have local rather than extensive spectral width; transfer layers have bounded local dimensions | Noncommuting bonds create splitting error; any trace monotonicity used must be proved for the exact powered product. Certification of its compressed powers is still required |
| Certify only the trace/purities with local or weighted residuals | Avoid charging every discarded component by the worst global operator norm | Requires a new rigorous contraction/error estimate. Small floating-point residuals alone do not certify it; cost depends on auxiliary rank and correlation length |
| Spatial-transfer contraction or a convergent expansion | Could control all lengths without evolving a 256-site operator | Need a checkable contraction/convergence condition at the target temperature. Assuming the desired gap or low-temperature clustering would be circular |
| Continue the present quadratic construction with symmetry alone | Existing checker and block arithmetic are available | Reduces arithmetic cost but leaves 28 squarings and potentially loose error amplification; not the first choice for another large run |

## Chosen diagnostic

Audit the scalar starting-point constraint and existing local-product experiments first.
Compare polynomial degree and squaring count with exact/outward-rounded bounds at N=256,
beta=64. Report the scalar improvement separately from the unresolved MPO construction and
compression costs. Use this to choose a meaningful-temperature small-chain benchmark before
another long-ring run. A route is not successful merely because its scalar remainder is small.

## Findings from the audit

### The tiny starting beta is a construction choice

The existing route centers K = H + N I and uses the elementary spectral enclosure
`-5N I <= K <= 5N I`. With `m = 2^power`, `delta = beta/(2m)` and polynomial degree d,
the implementation uses

```
x = 5N*delta
r >= exp(x)*x^(d+1)/(d+1)!
Q(K) = sum((-delta*K)^j/j!, j=0:d) + r I
log(trace inflation) <= 4m*r*exp(x).
```

The quadratic choice makes this restriction approximately proportional to
`(beta*N)^3/m^2` when x is small. Large N therefore forces many squarings even before
compression is considered. This is not a physical lower limit on the starting temperature.
Here the starting trace temperature is `2delta = beta/m`; distinguish it from delta.

A fresh 384-bit Arb acceptance check, using the existing `taylor_seed` generator, gave the
following first accepted powers at N=256, beta=64 and scalar inflation below 1001/1000:

| Degree | Subsequent squarings | Starting beta |
|---|---:|---:|
| 2 | 28 | 1/4194304 |
| 4 | 21 | 1/32768 |
| 8 | 17 | 1/2048 |
| 12 | 16 | 1/1024 |
| 16 | 15 | 1/512 |
| 32 | 13 | 1/128 |
| 64 | 12 | 1/64 |
| 128 | 11 | 1/32 |
| 256 | 10 | 1/16 |
| 512 | 9 | 1/8 |

Acceptance compared the rational log-inflation upper bound with an Arb lower enclosure of
`log(1001/1000)` and independently checked the remainder and log-inflation inequalities.
No Float64 comparison accepted a row. For the high degrees, a displayed Float64 inflation
can round to 1; it is not an assertion that the exact inflation is 1.

Higher-degree screening and the quartic SOS construction already existed; they are not new
methods discovered in this audit. The extended table is a negative check on a tempting
shortcut: merely raising the Taylor degree can exchange many squarings for a very expensive
initial MPO. No scalable degree-32 or degree-512 MPO certificate has been constructed.

### The local-gate route already escapes the tiny-beta regime

`src/ThermalBounds.jl` implements the powered even/odd product trace upper bound

```
Z_N(beta) <= Tr[(exp(-beta*A/m) exp(-beta*B/m))^m],   H = A+B.
```

This follows from Araki--Lieb--Thirring and the Lie product formula. It gives a rigorous
upper bound without a separate additive splitting-error allowance. Nevertheless, its
overestimate of the true Gibbs trace must be small enough to satisfy the seed target.
The existing four-site signed-residual/adjoint experiment starts at beta=1/8 and reaches
beta=4 in five doublings. Its compressed-backward certificate adds about 0.03011% over the
uncompressed Trotter trace (`results/adjoint-trace-compressed-n4.json`). It uses forward
rank 564 against full cut dimension 625, so this is not evidence of efficient scaling.

A fresh independent four-site calculation exposes an additional error budget. The spectrum
of H4 is exactly

```
E(u,v,s) = (s(s+1)-u(u+1)-v(v+1))/2, multiplicity 2s+1,
u,v = 0:4, s = abs(u-v):(u+v).
```

At beta=4 and shift 238061/50000, the exact shifted Gibbs trace is enclosed by an Arb ball
centered at 48.13451956631644. Fresh `trotter_upper` contractions give these ratios to it
(shown rounded to explanatory decimals):

| Trotter steps m | Local gate step beta/m | Trotter upper / true Gibbs trace |
|---|---|---:|
| 64 | 1/16 | 1.054501559780 |
| 128 | 1/32 | 1.013417771778 |
| 256 | 1/64 | 1.003341576520 |

Thus the earlier 0.03011% compression overhead sits on top of roughly 5.45% product-formula
overhead at m=64. Do not describe it as a 0.03011% approximation to the Gibbs trace. Both
budgets must be measured; the four-site ratios cannot be extrapolated into certified
256-site inequalities.

## Decision after the diagnostic

Do not launch another large quadratic-seed run merely because the contractions became
faster. Prioritize local-gate initialization plus a certificate for the final trace, using
the existing signed residual/adjoint machinery and charge-preserving arithmetic. First
separate splitting overhead, compression overhead, and verification cost on the four-site
and six-site references at beta of order 1 to 4. Then test dependence on length and retained
rank with the scalable MPO implementation. A high-degree polynomial is a comparison route,
not the default replacement.

The go/no-go question is whether the total upper bound stays useful as length increases;
neither a small scalar remainder nor a high-temperature prefix answers it. The old third-step
amplification near 161 is a bound from that representation and its norm estimates, not a
measured instability of the physical Gibbs operator. Its looseness needs its own diagnosis.

Full SU(2) reduction can reduce cost further, but does not by itself remove either source
of overestimation. A low-temperature expansion or spatial-transfer contraction would be a
larger conceptual advance only after a non-circular convergence/contraction condition is
established for this Hamiltonian.

## Reproduction and sources

Replay both tables and check all acceptance inequalities from the repository root:

```sh
OPENBLAS_NUM_THREADS=1 julia --startup-file=no --project=models/spin2 \
  models/spin2/scripts/polynomial/audit_thermal_start.jl /tmp/spin2-thermal-start-audit.json
```

Omit the output argument to print the audit without writing a file. The optional JSON records
strict rational enclosures, source and environment hashes, package versions, and explicit
false flags for compressed-MPO, large-ring thermal-seed, and uniform-gap certification.
The audit checks the ten reported squaring counts, independent 384-bit Arb remainder and
inflation inequalities, the full 625-state four-site multiplicity count, and decreasing
Trotter/Gibbs ratios above one for 64, 128, and 256 product steps. It uses exact integers and
rationals plus Nemo/Arb enclosures; no floating-point estimate accepts an inequality.
The run does not replay the historical compressed-MPO witnesses or modify `results/`.

Run with one BLAS thread and `--project=models/spin2`. The scalar diagnostic is reproducible
by evaluating `CertifiedGaps.CertifiedMPO.PolynomialBounds.taylor_seed(256, 64, p, d)` for
`p=0:32` and the table's degrees, checking its rational `log_inflation_upper` against
`log(Nemo.ArbField(384)(1001//1000))`. The remainder and log-inflation inequalities above
are the acceptance checks. The second table is obtained from
`Spin2Gap.ThermalBounds.trotter_upper(4, 4, m; shift=238061//50000)` divided by the exact
spectral sum above, both in Arb arithmetic. The code inspected was at commit
`22b2759`; only this memo is changed by the audit.

- [Existing scalar and SOS implementation](../../julia/src/CertifiedMPO/PolynomialBounds.jl).
- [Existing local-product reference](src/ThermalBounds.jl).
- [Existing trace correction benchmark](results/adjoint-trace-compressed-n4.json).
- [Audenaert, On the Araki--Lieb--Thirring inequality](https://arxiv.org/abs/math/0701129).
- [Al-Mohy and Higham, A New Scaling and Squaring Algorithm for the Matrix Exponential](https://eprints.maths.manchester.ac.uk/1300/).
  Overscaling is a known numerical issue, but that paper's norm-power improvement for
  nonnormal matrices is not a direct solution here: for Hermitian K in the physical metric,
  `||K^p||^(1/p) = ||K||`. The applicable lesson is to choose approximation order and scaling
  jointly, while retaining an independent certificate for every approximation.

## Outcome

The October 14 replay passed on Julia 1.12.4 with Nemo 0.57.0 and one BLAS thread.
The audit body took 0.43 seconds, excluding Julia startup, package precompilation, and
compilation before entering the function. Exact JSON parsing and all 42 recorded source
and environment hashes were checked separately. The new driver was formatted with
JuliaFormatter; `git diff --check` passed.

The audit identifies two distinct limitations: increasing polynomial degree reduces the
number of squarings without establishing affordable MPO construction, and the small
compression overhead at four sites excludes a larger Trotter overestimate. The selected
next research route remains local gates with final-trace certification and measured scaling.
The PR packages a reproducible diagnostic; a large-ring thermal seed and the uniform spin-2
gap certificate remain missing.
