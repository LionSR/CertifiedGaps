This draft adds a numerical BLBQ Haldane scaling survey and a constructive SU(2) highest-weight reduction under `certificates/blbq_haldane/`. It does **not** claim a uniform spectral-gap certificate.

Sweep checkpoint, 2026-10-11 07:04 UTC: 20/120 L14 sector files complete; S=2, k=4 active at 1m51s and 8.8 GiB/16 GiB. The preceding block completed in 7m53s. No timeouts or OOM events have occurred. Every sector has a 45-minute timeout. Remaining work is estimated at roughly 2–3 hours plus final audit, based on measured block timings; this is not a completion guarantee. Complete partition-function data still stops at L12.

- Exact Cartesian RP identities, including an independent realignment spectrum check, hold symbolically for the requested angle range.
- Complete dense spectra through L=12 cover all six angles. No eligible even seed passes; theta=0 is closest (n0=6, beta=3.25, balanced purity 0.328774).
- Numerical low-energy references at L=14 and 16 extend delta(n) through n=8. These truncated references are never used for partition functions. Their M=0 multiplet spacing equals the full gap only if the ground state is a unique singlet; this calculation does not certify that assumption.
- The new Casimir-polynomial highest-weight construction passes exact small projector checks and complete-spectrum comparisons through L=7. Large L=14 pilots succeeded. The full L=14 sector sweep is running as resumable, individually time-limited jobs; no incomplete sector set is treated as Z14.
- The theorem requires even n0, so full Z14 adds an odd-n0=7 diagnostic. The next eligible seed needs full Z16.

Validation and audit: 384 assertions passed in 93.8 seconds on current main. Julia 1.12.4; BLAS one thread. Current-main tests and namespace-migration checks are recorded separately from the preserved ThermalMPOGap snapshot. The core algorithms and noncore dependency entries are unchanged after the CertifiedGaps rename. Historical hashes are preserved; a pilot's end-of-run versus loaded-source hash correction is explicitly documented. Complete spectra have exact dimension/multiplicity totals and numerical trace-moment checks. No residual/orthogonality interval enclosure of eigenvalues has been produced, and no seed is certified. Extrapolations are conditional models; for theta <= -0.5 the required scale remains unresolved.

Reproduce with `OPENBLAS_NUM_THREADS=1 julia --project=julia certificates/blbq_haldane/test/runtests.jl`. The full commands, measured runtimes, scaling tables and limitations are in `REPORT.md`; migration details are in `MIGRATION.md`. `run_su2_grid.sh` caps every sector subprocess at 45 minutes and resumes completed sectors. `collect_su2.jl 14` refuses to assemble a partition function until all sectors exist.

Changes are confined to this certificate directory. No shared-core edits, CI, or merge are included. This remains a draft while the complete 14-site sweep and its final audit are pending.
