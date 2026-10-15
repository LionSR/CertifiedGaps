# Spin-1 large-D chain: uniform-gap certificates

H = Σⱼ Sⱼ·Sⱼ₊₁ + D Σⱼ (Sᶻⱼ)² on periodic rings of spin-1 sites, one site per cell, so all
lengths below count sites. Every conclusion is "unique ground state and Δ_L ≥ bound for every
even L ≥ n₀", under the repository's two-purity theorem. Bounds are rounded **down**; the exact
rationals are in the JSON files.

## Results

| D | Certified gap bound | Valid for | Seed (n₀, β₀) | Certificate |
|---|---:|---|---|---|
| 4 | **0.899143399** | every even L ≥ 6 | (12, 29/10) + finite bridge L = 6, 8, 10 | `certificate.json` |
| [3.2059998, 3.206] | **0.001001626** | every even L ≥ 6 | (6, 381/125), Weyl-transported | `certificate.json` |
| 4 | 0.694694192 | every even L ≥ 6 | (6, 9/4), ED seed | `results/certificate_D4_1.json` |
| 1603/500 = 3.206 | 0.002213456 | every even L ≥ 6 | (6, 381/125), ED seed | `results/certificate_D1603_500.json` |

Exact seed defects (εₛ ≥ 1 − Pₛ(n₀, β₀), εₜ ≥ 1 − P_τ(2n₀, β₀)):

| Seed | εₛ | εₜ |
|---|---|---|
| D = 4, n₀ = 12, β₀ = 29/10 | 8376656523711495 / 2^61 | 283096258190373 / 2^53 |
| D ∈ [3.2059998, 3.206], n₀ = 6, β₀ = 381/125 | 863115144109137 / 2^52 | 418622220052377 / 2^53 |
| D = 4, n₀ = 6, β₀ = 9/4 (stage 1) | 577405872923865 / 2^53 | 45632040552745 / 2^50 |
| D = 3.206, n₀ = 6, β₀ = 381/125 (stage 1) | 863088513139635 / 2^52 | 418496566764289 / 2^53 |

`certificate.json` also certifies, at D = 4, the larger seeds n₀ = 24, 48, 96 (bounds 0.660233530,
0.402885117, 0.063967113), the survey seed (6, 3) with bound **0.543812926** (the survey's
floating-point 0.543842 is neither used nor reproduced), and D ∈ [3.9979773, 4] with bound
0.001442884. For comparison, the finite-ring gaps are Δ₁₂ ≈ 2.26946 at D = 4 (Δ₈ ≈ 2.27392,
Δ₁₀ ≈ 2.27024) and Δ₁₂ ≈ 1.54533 at D = 3.206; the D = 4 bound is 0.396 Δ₁₂.

## Reflection positivity

Cartesian basis (x, y, z), Sₐ = −iGₐ with the real integer matrices

```
Gx = [0 0 0; 0 0 1; 0 -1 0]    Gy = [0 0 -1; 0 0 0; 1 0 0]    Gz = [0 1 0; -1 0 0; 0 0 0]
```

The checker verifies −h = Σₐ Gₐ ⊗ Gₐ (= Σₐ Gₐ ⊗ conj Gₐ) exactly in integer arithmetic, together
with the spin algebra and h = P − |Ω⟩⟨Ω|; the on-site term D(Sᶻ)² = D·diag(1, 1, 0) is real.

## How it works

**Stage 1 (exact diagonalisation, `certify.jl`, `seed.jl`).** Sᶻ, translation, spin inversion
at Sᶻ = 0 and conjugate-momentum pairing split 3¹² = 531441 states into 92 blocks (largest
5798). Every block is diagonalised densely, and each eigenvalue is enclosed by interval block
construction, rigorous residual and orthogonality norms, and a polar/Weyl bound. Interval log Z
at L = 6, 12 gives the exact rational ED seed n₀ = 6.

**Stage 2 (beyond ED, `transfer_certify.jl`).** From the same L = 6, 12 spectra:
transfer-operator squaring bounds Pₛ(12, β), and the exact coupling identity
P_τ(4n, β) = Pₛ(2n, 2β) P_τ(2n, β)² / Pₛ(2n, β)² bounds P_τ(24, β), **without building L = 24**.
This raises the D = 4 bound from 0.694694192 to 0.899143399 for even L ≥ 12; fresh verified
spectra at L = 8, 10 (and the stored L = 6 one) extend it to every even L ≥ 6. A centred Weyl
bound (0 ≤ Σ(Sᶻ)² ≤ L) transports a seed to an exact D interval. [METHOD.md](METHOD.md) proves
every inequality; [IDEAS.md](IDEAS.md) is the route memo.

## What is and is not rigorous

Rigorous (relative to the trusted base): the RP identity; sector completeness; every
eigenvalue enclosure at L = 6, 8, 10, 12; all log Z intervals; the defect bounds εₛ, εₜ (exact
rationals); transfer squaring, coupling and Weyl transport; and the recurrence bound from the
unmodified `CertifiedGaps.certified_gap_bound`. The trusted base is the repository's
two-purity theorem, Julia, IntervalArithmetic's outward rounding, and the standard IEEE
binary64 / BLAS dot-product error model used for the residual bounds (one BLAS thread).

Not rigorous, and never used as a premise: the floating-point scan below, the Δ₁₂ values used
as plot denominators (the certified markers do enclose the finite Δ₁₂, but it is not the
thermodynamic gap), every failure statement, and size extrapolations. D = 3.2059998 is the
smallest D certified by this search, not a threshold: the transport search fails at the next
10⁻⁷ grid point, and spatial lifts at D = 3.206 fail on β = 1.5:0.025:6 because thermal purity
is lost. These are failures of these estimates, not evidence of a vanishing gap.

The fast checker recomputes everything from the stored, hash-verified L = 6, 12 enclosures but
does **not** repeat their dense L = 12 eigensolves (35 min each); it does repeat the L = 8, 10
bridge eigensolves. A separate fresh dense replay of all L = 6, 12 blocks at both D reproduced
every eigenvalue centre and radius bit for bit, and an independent review of the stage-1
pipeline found no defect. The review is at
`git show 593299f:models/spin1_largeD/audit/INDEPENDENT_REVIEW.md`; the replay evidence was
removed before merging and its commit (`cf7bdd7`) is not on `main`. Re-running `certify.jl`
and comparing SHA-256 with `enclosure_hashes` is an equivalent replay.

## Numerical scan (floating point, not certified)

D = 1:0.25:4 and 3.205, 3.206; n₀ = 2…6 (odd n₀ only as diagnostics); β = 0.125:0.125:20,
refined to β = 2.75:0.002:3.75 at D = 3.205, 3.206. Summary in
[results/scan_summary.json](results/scan_summary.json); curve in
[results/transfer_curve.pdf](results/transfer_curve.pdf) (data with exact certified ratio
intervals in `results/transfer_curve.json`).

| D | best even seed (n₀, β₀) | numerical bound | Δ₁₂ | max_β min(Pₛ, P_τ), n = 6 | δ(6) = E₀(12) − 2E₀(6) |
|---:|---|---:|---:|---:|---:|
| 1 | none | — | 0.2554 | 0.335 | 0.3528 |
| 1.5 | none | — | 0.3418 | 0.375 | 0.3229 |
| 2 | none | — | 0.5840 | 0.494 | 0.2422 |
| 2.5 | none | — | 0.9452 | 0.655 | 0.1505 |
| 3 | none | — | 1.3641 | 0.799 | 0.0857 |
| 3.205 | none | — | 1.5444 | 0.839 | 0.0677 |
| 3.206 | (6, 3.048) | 0.00283 | 1.5453 | 0.839 | 0.0676 |
| 3.25 | (6, 3) | 0.0931 | 1.5845 | 0.847 | 0.0642 |
| 3.5 | (6, 2.625) | 0.3232 | 1.8095 | 0.888 | 0.0483 |
| 3.75 | (6, 2.375) | 0.5121 | 2.0381 | 0.916 | 0.0365 |
| 4 | (6, 2.25) | 0.6948 | 2.2695 | 0.938 | 0.0278 |

Missing seeds mean the tested grid failed, not that the gap closes. Dense L = 14 is out of
reach here: one Sᶻ = 1 momentum block has 41835 states, i.e. 28.0 GB for a single complex
matrix against a 16 GiB host ([results/feasibility.json](results/feasibility.json)); n₀ = 7 is
also odd and ineligible.

## Layout

| Path | Role |
|---|---|
| `model.jl` | model, RP identity, symmetry blocks, verified dense enclosures |
| `certify.jl`, `certify_parallel.jl` | stage 1: complete enclosed spectra (L = 12: ~35 min) |
| `seed.jl` | stage 1: ED seed certificates `results/certificate_D*.json` |
| `transfer_certify.jl` | stage 2: `certificate.json` and `results/finite_bridge.json` |
| `check.jl` | checker for both stages |
| `src/TransferCertificate.jl` | stage-2 proof arithmetic, bridge and checker |
| `src/spectra.jl` | byte-exact compact storage of the L = 12 spectra |
| `scripts/` | floating-point scan, β refinement, L = 14 feasibility, report/plot |
| `test/` | unit tests; `ed_seeds.jl` is the stage-1 replay used by `check.jl` |
| `results/` | spectra, certificates, bridge, scan summary, curve |

`model.jl`, `seed.jl`, `certify.jl` and `certify_parallel.jl` must stay byte-identical at the
top level: the stage-1 certificates and spectra record their SHA-256 under these paths (with
the core's former name `ThermalMPOGap`), and the checker verifies them. For that reason a few of
their lines exceed 100 columns. The L = 12 spectra are stored as `*.skeleton.json` plus
run-length-encoded `*.levels.bin`; `src/spectra.jl` rebuilds the original JSON bytes and the
checker compares their hash with the certificate.

## Reproduce

Julia 1.12.4; this directory is its own project depending on `../../julia`. From the
repository root:

```sh
export OPENBLAS_NUM_THREADS=1
julia --project=models/spin1_largeD -e 'using Pkg; Pkg.instantiate()'
julia --project=models/spin1_largeD models/spin1_largeD/check.jl            # ~20 s
julia --project=models/spin1_largeD models/spin1_largeD/test/runtests.jl    # ~25 s
julia --project=models/spin1_largeD models/spin1_largeD/transfer_certify.jl # ~15 s
```

Floating-point scan and plot (gitignored caches in `results/`; plotting needs `pdflatex`):

```sh
julia --project=models/spin1_largeD models/spin1_largeD/scripts/scan.jl              # ~4 min
julia --project=models/spin1_largeD models/spin1_largeD/scripts/scan.jl 1603//500 641//200
julia --project=models/spin1_largeD models/spin1_largeD/scripts/refine_beta.jl 1603//500 641//200
julia --project=models/spin1_largeD models/spin1_largeD/scripts/report.jl
```

Stage 1 is optional (the stored spectra are verified by hash). Each L = 12 run took
2117–2164 s with one BLAS thread; L = 6 takes under a second:

```sh
cd models/spin1_largeD
julia --project=. certify.jl 4 12                             # writes results/enclosed_*_L12.json
julia -t 4 --project=. certify_parallel.jl 1603//500 12
julia --project=. seed.jl 4 6 9//4                            # needs the full L = 6, 12 JSON
```

`results/feasibility.json` was produced on the original Linux host by the script now at
`scripts/feasibility.jl`; its recorded script hashes refer to the earlier path and content.
