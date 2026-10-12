# Independent review and fresh spectral replay

**Outcome: accepted within the stated trusted base** (repository two-purity theorem, Julia,
IntervalArithmetic outward rounding, the IEEE/BLAS dot-product error model).

`INDEPENDENT_REVIEW.md` is a separately delegated mathematical and code review of the
model-specific pipeline: exact RP and its normalisation, sector completeness, all three
spectral error sources (construction, residual, orthogonality), the polar/Weyl theorem, the
direction of the rational purity bounds, even-length interpolation and uniqueness. It found
no defect. `reviewer_checks.jl` adds an exact Cartesian-to-normalised-Sz bridge (129 checks)
and 256-bit residual/orthogonality diagnostics (90 checks); the diagnostics support, and do
not replace, the norm proof.

```sh
OPENBLAS_NUM_THREADS=1 julia --project=julia models/spin1_largeD/audit/reviewer_checks.jl
```

## Fresh dense replay (evidence removed from this directory)

The fast checker consumes stored eigenvalue radii. To close that gap, all L = 6 and L = 12
symmetry blocks for both D values were rebuilt from the exact model on the renamed core
(`CertifiedGaps` at `22b2759`), every dense eigensystem was recomputed, and new construction,
residual and orthogonality bounds were evaluated, without reading any stored spectrum. The
work was split into four single-BLAS-thread processes (charge groups {1,4,…,12} and {0,2,3},
the latter also redoing L = 6) taking 21.1–24.8 min each.

Every sorted eigenvalue centre and every radius was bitwise equal to the stored one, at both
lengths and both D, and both rational seed pairs and raw recurrence bounds reproduced exactly:

| D | β₀ | Reproduced bound | Displayed bound |
|---|---|---:|---:|
| 4 | 9/4 | 0.6946941929576096 | 0.694694192 |
| 1603/500 | 381/125 | 0.0022134568960786357 | 0.002213456 |

The per-block evidence (`audit/results/`, about 30 MB), the worker `replay_worker.jl`, the
assembler `assemble_replay.jl`, its checker and the logs were removed to keep the
certificate small; they are in the branch history at commit `cf7bdd7`. Re-running
`certify.jl` / `certify_parallel.jl` (see the top-level README) and comparing the output's
SHA-256 with `enclosure_hashes` in the certificates is an equivalent replay.

The fresh computation reuses the reviewed model and enclosure code; it is not a second,
independently written eigensolver, nor a formal verification of BLAS or of the theorem.
