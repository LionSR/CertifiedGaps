# Spin-one Heisenberg: Lean kernel feasibility of the seed

A feasibility study for checking the spin-1 Heisenberg seed inside the Lean kernel: a
Std-only Lake project with packed-`Nat` arithmetic (whole state vectors as one big integer),
`decide +kernel` benchmarks up to n = 12, kernel-checked Ritz/residual arithmetic witnesses
for H6, H8, H10, a Mathlib API audit and a review of the proposed spectral route. The packed
primitive is viable; the proposed end-to-end certificate and its 1.5–3-week estimate are not.

**No certificate.** No seed, eigenvalue enclosure or gap bound is claimed. The strongest
negative finding (Float64, non-certifying): with the untuned cut anchor "two open half chains
minus 4I", no rank 1 ≤ k < 3ⁿ satisfies the separation needed by the Weyl/counting argument
for H6 or H8, even with exact Ritz vectors; the best margin is −1.383 (n = 6) and −1.445
(n = 8).

| n | Lean check | wall time (M5) |
|---|---|---:|
| 6 | dense 128-bit matvec + dot (Bench6) | 0.5 s |
| 8 | dense 128-bit matvec + dot (Bench8) | 0.5 s |
| 10 | dense 128-bit matvec + dot (Bench10) | 1.3 s |
| 12 | dense 128-bit matvec + dot (Bench12) | 26.4 s |
| 6, 8, 10 | Ritz ≤ −8.495447, −11.179504, −13.548710 (arithmetic witnesses) | 8.5 s (whole build) |

Times are wall-clock per fresh Lean process including imports (the imports-only baseline is
0.4–3.2 s), Apple M5, Lean v4.35.0-rc3.

## Model

Spin-1 antiferromagnetic Heisenberg ring H = Σⱼ Sⱼ·Sⱼ₊₁ on n sites, in the standard Sᶻ basis
(base-3 digits sⱼ = mⱼ + 1, index Σ sⱼ 3ʲ), so that exchange hops are integers. `hApply`
supports a π twist about z on the closing bond (sign flip of its transverse hops). The power
filter uses nI − H_M with H_M the Marshall conjugate of the untwisted even ring (nonnegative
entries); the Ritz witnesses use the trial v = (nI − H_M)^q e_Néel with
(n, q) = (6, 24), (8, 32), (10, 16). Lengths count sites.

## Reflection positivity

Not re-proved here. The certified spin-1 route (`models/spin1/`, `lean/`) uses the Cartesian
basis. The unitary bridge from that Hamiltonian to this Sᶻ-basis integer matrix, and the
twisted-ring extension, are open; METHOD.md §4 (Step 4) states the required lemmas.

## Result

Kernel-checked (`decide +kernel`, Std only, standard axioms):

| n | degree q | lane bits | Ritz witness θ ≤ | residual witness ≤ |
|---|---:|---:|---:|---:|
| 6 | 24 | 262 | −8495447/10⁶ | 274967/10⁶ |
| 8 | 32 | 345 | −11179504/10⁶ | 279904/10⁶ |
| 10 | 16 | 222 | −13548710/10⁶ | 632355/10⁶ |

These are integer inequalities `ritzCheck` on the moments G = vᵀv, A = vᵀFv, C = (Fv)ᵀFv
computed by the packed arithmetic. They become physical eigenvalue statements only after
the open packing-correctness obligations (`Obligations.lean`: `DotCorrect`, `HApplyCorrect`,
`FilterCorrect`, ...) and the Cartesian bridge are proved; the residuals are far too loose
for the thermal tables anyway. Also proved: elementary radix-packing lemmas
(`PackingProofs.lean`) and a conditional one-Ritz-value counting step (`Counting.lean`).

Findings (details, tables and the 8.5–15k-line implementation budget in METHOD.md):

- **Anchor.** The untuned cut comparison fails the separation test at n = 6 and 8 by more
  than one energy unit (`results/anchor_separation_noncertifying.tsv`). A tuned comparison
  with a positive counted margin must be found before a recursive Lean implementation is
  worth budgeting.
- **Compact packing.** The Sz = 0 sector at n = 12 has 73789 states but 3804 distinct
  rank shifts (3744 from the closing bond), refuting the "about 50 masked shifts" plan for
  straightforward compact packing (`results/structure.tsv`). Full packing (531441 lanes) is
  7.2x larger but has fixed shifts.
- **Cost.** At about 30 s per n = 12 matvec (Linux run), k = 50 vectors and q = 200 filter
  steps would take about 83 h for one n/twist before Gram work (the M5 re-run is ~3x faster).
- **Mathlib (pinned c55e6e78).** Present: Hermitian spectral theorem, descending ordered
  eigenvalues, trace = Σ eigenvalues, extremal Rayleigh eigenvalues, PSD/Loewner order,
  Jensen and `convexOn_exp`. Missing: Courant–Fischer min–max, Weyl inequalities for
  ordered eigenvalues, Cauchy interlacing/indexed Rayleigh–Ritz. Estimated 2.8–4.9k Lean
  lines for a focused spectral foundation.

## Reproduce

From this directory, with `lake` (elan) and Julia 1.12 on PATH, `OPENBLAS_NUM_THREADS=1`.
Times on an Apple M5.

```sh
julia --project=. -e 'using Pkg; Pkg.instantiate()'
(cd src && lake build)                               # 8.5 s from clean; no Mathlib needed
julia --project=. test/runtests.jl                   # 2 s, exact generator tests
julia --project=. check.jl                           # ~2 min: hashes, regeneration, clean build
julia --project=. check.jl --full                    # ~3 min: + axiom audit, every Bench file
julia --project=. scripts/generate.jl                # 105 s: H{6,8,10}.lean + exact_witnesses.tsv
julia --project=. scripts/witnesses.jl               # results/witnesses.json (hashes)
julia --project=. scripts/experiments.jl             # 5 s: Float64 structure/anchor tables
julia --startup-file=no scripts/benchmark.jl 300     # ~2.5 min: results/benchmarks.tsv
scripts/audit_mathlib.sh                             # needs lean/.lake (or MATHLIB_DIR)
scripts/fetch_references.sh                          # paper excerpts into references/
```

The project builds against Std alone. `audit_mathlib.sh` only reads Mathlib source; by
default it uses the checkout of the repository's `lean/` project, which pins the same
toolchain and Mathlib commit. Raw Lean benchmark output goes to the ignored `results/bench/`.

## What is rigorous

- Trusted base: the Lean 4 kernel (v4.35.0-rc3) evaluating `decide +kernel` on `Nat`
  (GMP-accelerated kernel arithmetic); no `native_decide`, no `sorry`, no new axioms.
  `#print axioms` on all 11 audited declarations gives subsets of `propext`,
  `Classical.choice`, `Quot.sound` (`src/AxiomAudit.lean`, `results/axioms.txt`).
- Kernel-checked: the integer statements in the table above, the benchmark identities,
  the packing lemmas, the conditional counting step.
- Not proved: correctness of packed masks, shifts, reversal, dot and H application against
  a matrix (stated as `Prop` definitions in `Obligations.lean`, never assumed); the
  Cartesian ↔ Sᶻ bridge; any ordered-eigenvalue, Weyl or counting theorem for H; the anchor.
- Floating point only: `results/*_noncertifying.tsv` and `results/structure.tsv` counts
  (exact enumeration, but in Julia, not Lean).
- The Julia generator is untrusted; its agreement with independent sparse rows is a
  regression check only.

## Layout

| path | content |
|---|---|
| `src/` | Lake project (Std only): `lakefile.toml`, `lean-toolchain`, `lake-manifest.json` |
| `src/SeedPrototype/` | library: `Packed`, `Optimized`, `PackingProofs`, `Obligations`, `Counting`, `Tests`, generated `H6/H8/H10` |
| `src/Bench/` | `decide +kernel` benchmark files (not in the default target) |
| `src/AxiomAudit.lean` | `#print axioms` for the main declarations |
| `scripts/` | `generate.jl`, `witnesses.jl`, `experiments.jl`, `benchmark.jl`, `audit_mathlib.sh`, `fetch_references.sh` |
| `check.jl`, `test/runtests.jl` | checker and fast tests |
| `results/witnesses.json` | exact witness inventory with source SHA-256 (not a certificate) |
| `results/exact_witnesses.tsv` | witnesses with the exact moments G, A, C |
| `results/benchmarks.tsv`, `benchmarks_original_linux.tsv` | kernel timings (M5 re-run, original) |
| `results/structure.tsv`, `*_noncertifying.tsv` | packing structure and Float64 anchor/level scans |
| `results/axioms.txt`, `mathlib-audit.txt`, `mathlib-negative-search.txt` | audits |
| `results/paper-source-sha256.txt` | hashes of the cited paper excerpts |
| `METHOD.md` | Mathlib audit, trust boundary, mathematical review, benchmarks, budget |

## History

PR #21 (prototype and audit; written where the Mathlib cache was unreachable, hence Std).
Sibling studies of the same seed live in `../spectrum` (seed spectra, precision, kernel cost
model), `../anchor` (local PSD comparisons for the anchor problem found here) and
`../seed_design`. Cleanup: standard layout with the Lake project under `src/`; dropped the
unused Mathlib requirement; removed build logs, the `cache-get` log and Linux-only tooling;
folded `REPORT.md` into this README and METHOD.md; re-ran every check and benchmark on an
Apple M5.
