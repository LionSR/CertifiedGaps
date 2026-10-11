Adds complete finite-spectrum gap certificates for the spin-1 antiferromagnet with single-ion anisotropy, H=Σ Sⱼ·Sⱼ₊₁+DΣ(Sᶻⱼ)². The smallest successful sampled D is 3.206; the neighboring D=3.205 sample has no seed on the refined temperature grid. This is a search threshold for n₀≤6, not a proven minimum over real D or a phase boundary.

| D | Seed (n₀, β₀) | Certified uniform lower bound, rounded down |
|---|---|---:|
| 4 | (6, 9/4) | 0.694694192 |
| 1603/500 = 3.206 | (6, 381/125) | 0.002213456 |

Both bounds apply to **every even L≥6 one-site cells, hence sites**, under the repository's two-purity theorem. See [REPORT.md](https://github.com/LionSR/haldane-gap-diagrams/blob/codex/spin1-large-d-certificate/certificates/spin1_largeD/REPORT.md) for the scan, ED reference gaps, bound/gap curve, runtimes, negative-grid diagnostics and extrapolations; [METHOD.md](https://github.com/LionSR/haldane-gap-diagrams/blob/codex/spin1-large-d-certificate/certificates/spin1_largeD/METHOD.md) gives the enclosure proof.

The exact integer Cartesian RP identity is checked independently against P−|Ω⟩⟨Ω|. Certification uses every dense symmetry block at L=6 and L=12, with multiplicities summing to 3^L (531441 levels at L=12). Interval matrix construction plus residual and orthogonality norm bounds feed a polar-decomposition/Weyl eigenvalue enclosure. Outward partition-function intervals give rational defect bounds stored as strings and passed to the unmodified core recurrence. The two large runs took 36.1 and 35.3 minutes; no Lanczos spectrum enters either certificate.

Only `certificates/spin1_largeD/` changes. No core, other model, CI or Lean changes. The explicit requested certificate location is retained despite the upstream directory migration.

## Validation and provenance

- Original-core certificate checker: **123/123**, 13.5 s wall time.
- Current `CertifiedGaps` checker: **134/134**, 19.6 s including precompilation; both bounds reproduced exactly.
- Pinned historical replay: **123/123**, 18.9 s wall time.
- Independent full dense comparisons through L=6 and enclosure checks through L=4: **20/20**, repeated on the current core.
- Exact translation/spin-inversion/charge checks and orbit counts through L=12: **3278/3278**, repeated on the current core.

The dense jobs used `ef0a3e6cf54ec56dcef5e8a37963d0184f6fd31b` (`ThermalMPOGap`). Main subsequently renamed the package at `22b2759e0790197b83e2571b465aa8dc69ebc0be`. Original source bytes are archived under `provenance/legacy/`; historical hashes were not rewritten. Imports were migrated only inside the certificate directory. The new checker verifies package-name-only changes to model/drivers, immutable historical core objects, unchanged recurrence source, and the stored spectra/seed replay. `results/validation-current.json` records the separate validation.

## Audit limitations

This remains a draft for scientific review. The eigenvalue proof assumes IEEE binary64 arithmetic and the documented standard BLAS dot-product error model; it explicitly bounds BLAS rounding rather than treating floating eigensolver output as exact. The fast checker verifies provenance, stored spectral-enclosure accounting, partition intervals and the rational recurrence; it does **not** reconstruct every eigenvector/residual. Full dense replay is needed to repeat that stage. The upstream rename validation is not a second L=12 ED run or a blanket audit of the entire core.

The search uses numerical partial spectra only for screening. Unsuccessful grid points, ED reference gaps and ratios, interpolation in the plot, and size extrapolations are **numerical, not certified**. No nonexistence result over all D or β is claimed. L=14 dense ED was not feasible on this 16 GiB host: a charge-one momentum block alone needs 28.0 GB for one complex dense matrix. The supplied D=4, β=3 survey value differs from the numerical rerun by about 2.17e−5; it was not used as certificate input.

## Reproduction

Julia 1.12.4, instantiated `julia/` dependencies, from the repository root:

```sh
export OPENBLAS_NUM_THREADS=1
julia --project=julia certificates/spin1_largeD/test/runtests.jl
julia --project=julia certificates/spin1_largeD/test/small.jl
julia --project=julia certificates/spin1_largeD/test/symmetry.jl
julia --project=julia certificates/spin1_largeD/check_pinned.jl
```

The pinned replay needs the original core commit in Git history and installed dependencies; it uses a temporary workspace and changes no checkout files. For full dense reproduction, use a disposable checkout to preserve the recorded results:

```sh
julia --project=julia certificates/spin1_largeD/certify.jl 4 6
julia --project=julia certificates/spin1_largeD/certify.jl 4 12
julia --project=julia certificates/spin1_largeD/seed.jl 4 6 9//4
julia --project=julia certificates/spin1_largeD/certify.jl 1603//500 6
julia --threads=2 --project=julia certificates/spin1_largeD/certify_parallel.jl 1603//500 12
julia --project=julia certificates/spin1_largeD/seed.jl 1603//500 6 381//125
```

Historical full replay instructions, plotting commands, source hashes and the preserved original prompt are in the directory README and results.
