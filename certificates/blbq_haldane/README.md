# BLBQ Haldane finite-size survey

See [REPORT.md](REPORT.md) for results, exact Cartesian reflection positivity, numerical
purity and energy-mismatch tables, resource limits, and conditional seed-scale forecasts.

This directory is a **numerical negative grid survey**, not a uniform-gap certificate.
Complete dense spectra reach 12 sites; low-energy numerical references extend to 16 sites.
The numerical SU(2) highest-weight construction is implemented and validated against
small complete spectra. Large 14-site pilots succeeded; the complete sector sweep is
being collected. Selected sectors must not be used as a complete partition function.

From the repository root, with Julia 1.12.4 and the core dependencies installed:

```sh
OPENBLAS_NUM_THREADS=1 julia --project=julia certificates/blbq_haldane/test/runtests.jl
```

`run.sh` additionally selects the locally downloaded runtime and private depot if available.
All results JSON uses `CertifiedGaps.Results.write_result`; exact rational tokens are strings.
The failed numerical seed checker does not validate numerical defects as rigorous enclosures.
