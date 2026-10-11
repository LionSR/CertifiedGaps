# BLBQ Haldane finite-size survey

See [REPORT.md](REPORT.md) for results, exact Cartesian reflection positivity, numerical
purity and energy-mismatch tables, resource limits, and conditional seed-scale forecasts.

This directory is a **numerical negative grid survey**, not a uniform-gap certificate.
Complete dense spectra reach 14 sites for all six angles; low-energy numerical references
extend to 16 sites. The SU(2) highest-weight construction is validated against small
complete spectra. All 120 L14 sector entries are present, with exact integer multiplicity
counts totaling 3^14 states. L14 adds the odd-n0=7 diagnostic; it does not supply an
admissible even seed. The next eligible seed requires complete L16 data, beyond this
implementation's memory budget.

From the repository root, with Julia 1.12.4 and the core dependencies installed:

```sh
OPENBLAS_NUM_THREADS=1 julia --project=julia certificates/blbq_haldane/test/runtests.jl
```

`run.sh` additionally selects the locally downloaded runtime and private depot if available.
All results JSON uses `CertifiedGaps.Results.write_result`; exact rational tokens are strings.
The failed numerical seed checker does not validate numerical defects as rigorous enclosures.
