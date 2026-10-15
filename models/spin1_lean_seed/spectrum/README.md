# Spin-1 seed spectra, filtering precision and Lean kernel costs

**No certificate; non-certifying numerical audit.** This directory measures what a Lean-checked
seed for the spin-1 Heisenberg certificate (`models/spin1/`) would have to compute: how many low
states each twisted ring keeps, how accurately they must be known, whether integer Chebyshev
filtering reaches that accuracy, and what the arithmetic would cost in the Lean kernel. All
spectra are Float64; the only kernel-checked statements are arithmetic benchmarks.

| source | β | rings | states kept at n_max (I / P / C) | adaptive Q48 Hv | Lean vector + Gram forecast |
|---|---:|---|---|---:|---|
| 1 | 49/4 | n = 4–12, I/P | 26 / 34 / – at n = 12 | 34,388 | 38–83 min (plan 1–4 h) |
| 2 | 21/2 | n = 6, 8, 10, I/P/C | 22 / 30 / 26 at n = 10 | 17,800 | 5–11 min (plan 10–30 min) |

Two premises of the proposed route fail (direct compressed-Sz shifts; the naive spectral
anchor); see "Result".

## Model

Spin-1 Heisenberg ring `H_n = Σ_j S_j·S_{j+1}` in the Sz product basis, one site per cell, with
the closing bond (n, 1) twisted by a rotation about z of 0 (I), π (P) or 2π/3 (C): its flip
terms carry 1, −1 or exp(±2πi/3). Bond eigenvalues are −2, −1, 1, so n bounds the spectrum above.
Traces use the shift a = 700741/500000 of `models/spin1/`, `Ẑ = Tr exp(−β(H + a n))`.
Source 1 is β = 49/4 (the original (36, 49/4) initialization) on n = 4–12 with I/P; source 2 is
β = 21/2 on n = 6, 8, 10 with I/P/C. Each Sz sector is a separate matrix; negative Sz sectors
are inferred from spin inversion (for C combined with complex conjugation).

## Reflection positivity

Not re-examined. The existing spin-1 proof uses the Cartesian RP basis (`models/spin1/`); this
directory only computes spectra in the Sz basis, which needs no RP structure.

## Result

All numbers are floating point (`METHOD.md` has the derivations and full tables).

**Low-state counts.** k(c) counts eigenvalues below E₀ + c over all Sz sectors with
multiplicity; c* is the smallest cutoff with dimension tail `(3ⁿ − k(c)) e^{−βc} / Q < 10⁻⁶`
(Q the ground-relative partition function). At the largest rings:

| n | twist | E₀ | k(1) | k(1.5) | k(2) | k(2.5) | k(3) | c* / k*, β = 12.25 | c* / k*, β = 10.5 |
|---|---|---:|---:|---:|---:|---:|---:|---|---|
| 10 | I | −14.0941 | 4 | 9 | 21 | 26 | 67 | 2.0242 / 21 | 2.36088 / 22 |
| 10 | P | −13.9821 | 6 | 11 | 20 | 36 | 74 | 2.02395 / 20 | 2.36048 / 30 |
| 10 | C | −14.0138 | 4 | 9 | 20 | 35 | 70 | 2.02412 / 20 | 2.36076 / 26 |
| 12 | I | −16.8696 | 4 | 15 | 22 | 48 | 91 | 2.20334 / 26 | 2.56956 / 58 |
| 12 | P | −16.7976 | 6 | 14 | 24 | 56 | 99 | 2.20335 / 34 | 2.56974 / 60 |

"About E₀ + 2" understates the cutoff at n = 12: the c = 2 tail is 1.21·10⁻⁵ (β = 12.25) and
3.96·10⁻⁴ (β = 10.5). Per-sector counts are in `results/sectors.csv`; every ring is in
`METHOD.md`. Precision is not monotone in n: an absolute trace error 3·10⁻⁵ needs eigenvalue
accuracy 1.96·10⁻⁸ at n = 4/I (Ẑ = 124.885) but only 1.29·10⁻⁶ at n = 12/I
(`results/precision.csv`).

**Filtering precision.** Block Chebyshev filtering (block k + 8, QR and Rayleigh–Ritz after each
cycle) from computational-basis states, with the recurrence in exact Int128 fixed point (Q48 =
64-bit lanes, Q112 = 128-bit) and Float64 QR. At n = 12/I/Sz = 0: degree-32 Q48 plateaued at
residual **2.99·10⁻⁶** after 40 cycles; Q112 reached **2.34·10⁻⁹** in 4 cycles; degree-16 Q48
reached **3.29·10⁻⁸** in 8 cycles. Uniform degree-32 Q48 misses the linear target in two sectors
per source (Sz = 0 of n = 12/I, P and of n = 10/I, P); Q112 degree 32 passes everywhere.
Degree-64 Q48 stalls (n = 12, Sz = 0: 0.0666 for P, 0.269 for I). With the cheapest successful
measured degree per sector, every active sector meets both residual targets: source 1 44/44
(Q48 linear 34,388 Hv, quadratic 26,562), source 2 24/24 (17,800 and 12,486). Final Gram cost
is 1,977 (source 1) and 2,100 (source 2) scalar products.

**Packing.** Direct masked shifts on the compressed Sz basis need many displacement masks, not
about 50 (`results/packing.csv`, both directions per bond, diagonal extra):

| n | Sz | dimension | full lanes 3ⁿ | compressed masks | full-basis masks |
|---|---:|---:|---:|---:|---:|
| 8 | 0 | 1,107 | 6,561 | 126 | 16 |
| 10 | 0 | 8,953 | 59,049 | 624 | 20 |
| 12 | 0 | 73,789 | 531,441 | 3,816 | 24 |
| 12 | 1 | 69,576 | 531,441 | 9,884 | 24 |
| 12 | 2 | 58,278 | 531,441 | 13,620 | 24 |

**Cost model** (`results/cost_totals.csv`, `results/runtime_estimates.csv`; Q48, linear target,
adaptive degrees; vector generation, modelled mixing and final Gram only):

| source | layout | Hv | bit-pass volume | measured-rate model (min) |
|---|---|---:|---:|---|
| 1 | full ternary | 34,388 | 4.15·10¹³ | 41.3–82.5 |
| 1 | full Hv + compacted Gram | 34,388 | 4.30·10¹³ | 38.1–74.7 |
| 1 | compressed direct | 34,388 | 9.58·10¹⁴ | 959–1758 |
| 2 | full ternary | 17,800 | 4.64·10¹² | 5.73–10.6 |
| 2 | full Hv + compacted Gram | 17,800 | 4.82·10¹² | 5.26–9.39 |
| 2 | compressed direct | 17,800 | 2.43·10¹³ | 49.5–69.9 |

Source 2 needs about 9× less full-layout work than source 1, so lowering n_max is the main lever.
The rates come from Lean 4.35.0-rc3 `decide +kernel` benchmarks on a Xeon Platinum 8573C: one
4M-bit product/division 0.91 s with startup, about 0.016–0.024 s per 4M-bit multiplication in
bulk, 5.6–9.4·10⁻¹¹ s per mask/shift/XOR bit pass, about 31k small-loop iterations/s.

**Anchor (the open problem).** Cutting two bonds (`h ≥ −2`) gives
`H ≥ H_open(⌊n/2⌋) + H_open(⌈n/2⌉) − 4`. At n = 12/I/Sz = 0 the anchor ρ = −14.6662 must
exclude all but k = 10 levels; the comparison operator has **242** levels below ρ and its 11th
level is −17.0810 (deficit 2.4148). Over all sectors it has **968** levels below ρ against 26
kept, and it cannot exclude the k = 0 sectors Sz = 3, 4. At n = 10/C/Sz = 0 (source 2): k = 8
against 165 (deficit 2.5824). Tuned boundary terms or a homotopy were not tested.

## Reproduce

Julia 1.12, one BLAS thread, from this directory (`export OPENBLAS_NUM_THREADS=1`):

```sh
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. test/runtests.jl           # 1,459 checks, about 20 s
julia --project=. check.jl                   # the 450 data checks alone, about 12 s
```

Regeneration (outputs go to `results/`, or to `$SPIN1_SPECTRUM_RESULTS` if set; spectrum and
filter JSONs act as a cache, so move them aside to recompute). Times are on an Apple M-series
laptop unless marked:

```sh
julia --project=. scripts/spectra.jl            # all 21 rings, 2.6 min
julia --project=. scripts/analyze.jl            # sectors/traces/precision.csv, eigenvalues.csv
julia --project=. scripts/filters.jl float      # Float64 degrees 16/32/64 (not re-timed)
julia --project=. scripts/filters.jl fixed      # Q48/Q112 degree 32 (n = 10 alone: 34 s)
julia --project=. scripts/filters.jl sweep      # Q48/Q112 degrees 16/64 at n = 10, 12
julia --project=. scripts/packing.jl            # packing.csv/json
julia --project=. scripts/anchor_probe.jl       # anchor_probe.csv
julia --project=. scripts/summarize.jl          # filter_runs/curves, cost_*.csv, summary.json
julia --project=. scripts/runtime_estimates.jl  # runtime_estimates.csv
julia --project=. scripts/provenance.jl         # provenance.json (hashes)
```

`analyze`, `packing`, `anchor_probe`, `summarize` and `runtime_estimates` together take about
10 s. The original runs took minutes per stage on a Xeon Platinum 8573C (5 CPUs, independent
processes). `summarize.jl` reads the raw `filter_s*.json` if present; otherwise it rebuilds the
cost tables from the runs stored in `summary.json` (both paths give byte-identical `cost_*.csv`).
Three stored degree-32 Q48 failures (n = 12/I and n = 10/I, P at Sz = 0) ran 40 cycles; the
current default stops at 16 or a detected plateau (8 cycles for the two n = 10 runs on
regeneration, with the same failure). Regenerated spectra agree with the stored ones to
2.5·10⁻¹³ with identical counts; a few near-threshold cycle counts may differ across BLAS builds.

The optional Lean benchmarks need Lean 4.35.0-rc3 at `runtime/lean-4.35.0-rc3-linux/bin/lean`
(ignored): `scripts/kernel_bench.jl`, `scripts/kernel_bulk.jl`, and `scripts/KernelBench.lean`.
They write temporary Lean files under `runtime/` only; no `native_decide`, `sorry` or extra axiom
(`#print axioms` reports none or `[propext]`; in `KernelBench.lean` it names the definition
`packedA`, not the anonymous `example`, so that line does not audit the checked statement).

## What is rigorous

Nothing here bounds a partition function. Trusted: the Int128 floor-multiplication primitive
(tested against BigInt) and the kernel-checked Lean arithmetic benchmarks. Floating point only:
all spectra (dense LAPACK through n = 8; ARPACK above, with residuals < 10⁻⁸ and guard levels
beyond E₀ + 3.15), counts, cutoffs, filter residuals, QR, anchor probes and cost forecasts. A
small residual does not show that Lanczos missed no eigenvalue, so every count and anchor is
conditional on the observed ordering; the Ritz-value comparison in `check.jl` uses the same
ARPACK list. The trace check against `models/spin1/results/source{1,2}.json` (27 centres within
3·10⁻⁵) is a cross-check only. Regenerated CSVs can differ from the stored ones in the last digit
on other platforms (libm `exp`, LAPACK).

## Layout

- `src/common.jl`: basis, sparse twisted Hamiltonian, ARPACK/dense low spectrum, I/O.
- `src/filtering.jl`: Chebyshev coefficients, Float64 and fixed-point filters, sector runs.
- `scripts/`: the drivers listed under "Reproduce", plus the Lean benchmark files.
- `check.jl`: data consistency checks (not a certificate checker); `test/runtests.jl` runs the
  integer-primitive tests and then `check.jl`.
- `results/spectrum_<n>_<twist>.json`: per-sector Float64 spectra (complete through n = 8).
- `results/analysis.json`, `summary.json`: inputs of `check.jl` and of the filter/cost scripts.
- `results/*.csv`: the tables above (sectors, traces, precision, filter runs and curves, packing,
  cost totals, runtime estimates, anchor probe, kernel timings, reference trace checks);
  `kernel-*.txt` are the Lean benchmark logs; `provenance.json` hashes sources and data.
- Not in git (regenerable, `.gitignore`): `results/filter_s*.json` (424 raw filter runs, 3.6 MB),
  `eigenvalues.csv` (1.6 MB, from `analyze.jl`), `cost_sectors.csv` (0.7 MB, from
  `summarize.jl`), `dense_validation.csv` (written by `check.jl`).

The `source_hashes` inside each `spectrum_*.json`, `analysis.json` and `summary.json` record the
flat-layout scripts at generation time (some are development versions); `provenance.json` hashes
the current layout and keeps the original generation machine's metadata.

## History

PR #22 (`codex/spin1-seed-spectrum`), after the route memo in `IDEAS.md`. Cleanup: moved to the
standard layout (`data/` → `results/`, scripts into `src/` and `scripts/`, `REPORT.md` →
`METHOD.md`, `report.jl` dropped since README and METHOD are now maintained by hand,
`validate.jl` merged into `check.jl`), compacted the JSON (parsed contents identical), and
removed the regenerable dumps (tracked data 7.8 MB → 1.8 MB).
Follow-ups: `../anchor/` replaces the naive two-cut comparison here by exact local PSD
comparisons (its "bare" count, 191 at n = 12/I/Sz = 0, uses its own cutoff ρ ≈ −14.891, not the
ρ = −14.666 above), and `../seed_design/` reduces the finite seed inputs.
