# Validation record

Lean `4.35.0-rc3` (commit `470d5ce1400764999581fd26d5d72b00d990b0f4`),
Mathlib `c55e6e786f49471c72fbddbec5415808896aec1e`, Julia `1.12.1`, BLAS threads `1`.

| Check | Result | Evidence |
|---|---|---|
| Packed integer checks, n=6 and n=12 | PASS | `check6`, `check12`, using `decide +kernel`; axioms `[propext]` |
| Direct scalar integer check, n=6 | PASS | `direct6`; axioms `[propext]` |
| Direct scalar integer check, n=12, 64 blocks | PASS | `direct12_chunked`; axioms `[propext]`; module built in 123 s |
| Coverage of all 4096 binary states | PASS | `coverage12`; no axioms |
| Classical diagonal counts, n=6 and n=12 | PASS | `anchor6`, `anchor12`; no axioms |
| Exact rational seed endpoint arithmetic | PASS | `spatial_arithmetic`, `thermal_arithmetic`; standard three axioms |
| Weighted residual existence lemma | PASS | `exists_near_of_weighted_residual`; standard three axioms |
| False eigenvalue/vector rejection | PASS | Two kernel-checked negative tests; axioms `[propext]` |
| Julia exact replay and mutation tests | PASS | 8/8 tests; about 9 s after startup; diagnostic only |
| All-model inventory | PASS | 43 rows, dimension assertions and grid-point identity assertions |
| Numerical two-level feasibility | PASS, non-certifying | Both defects pass; all four proposed log bounds fit the source endpoints |
| Integrated Lake build and exact RPRing bridge | PASS | Final successful build: 2727 jobs, 5.6 s with dependencies already built |
| Literal first instance of the original 30-point seed statement | PASS, conditional | `first_point_seed_instance`; only the standard three axioms |
| Final axiom and source audit | PASS | All 17 declarations; protected Lean source hashes and both pins match |

The compact packed module built in about 53 s during a concurrent Mathlib source build.
Profiling separately showed approximately 0.3 s for the 12-site packed arithmetic, 6 s for
compact literal elaboration, and 52 s for direct classical enumeration. These are observed
times on this environment, not portable performance guarantees.

The single-expression direct 12-site attempt was deliberately stopped after a snapshot showed
7.4 GB RSS and 38 CPU seconds without completion. Separating the scalar work into independently
checked blocks produced a successful check. This failure is retained as a cost observation,
not suppressed or counted as a successful test.

`validate.jl` checks the final build, all 17 printed axiom lists, Julia tests, the exact Lean and
Mathlib pins, and hashes of the protected Lean sources. It rejects `sorryAx`, native reduction
axioms, and any axiom outside `propext`, `Classical.choice`, `Quot.sound`. It writes
`validation.json`; logs stay in the ignored `logs/` directory.
The final axiom output is also retained in [AXIOMS.txt](AXIOMS.txt). Final build timing is
incremental; the 53 s and 123 s module measurements above were actual initial proof checks.

**The seed conclusions remain conditional on `LogEnclosures`.** Passing an axiom audit does
not prove those four hypotheses. The matrix identification, norm interpretation, spectral
counting, and scalar exponential-enclosure proof are not claimed as complete.
