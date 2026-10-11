# Preserved pre-migration evidence

The original complete spectra through L=12 and low-energy references through L=16 were
computed with `ThermalMPOGap` on repository commit `ef0a3e6` (Julia 1.12.4). Original result
metadata and source hashes remain in those result files. The L=12 JSON is losslessly
compressed; `results/archive.json` records both its original and compressed hashes.

`thermal_validation.json`, `thermal_scaling.json`, `thermal_environment.json`, and
`thermal_rp.json` preserve the first completed survey and its validation unchanged.
`thermal_source/` preserves the original model, sector solver, scan, validation driver and
test file; every saved file was checked against its original recorded SHA-256 before saving.
For old-snapshot reproduction, use the core at ef0a3e6 and restore these files to the
certificate root; decompress the L=12 JSON to its original filename for the old checker.
The compressed archive hashes to the original JSON, rather than being a recomputation.

`su2_householder.jl` preserves the first numerical highest-weight implementation exactly.
Its SHA-256 matches the L=12 pilot's recorded source hash. It is a source snapshot; its
relative includes assume placement at the certificate root. The later CholeskyQR2
implementation is separately validated and benchmarked.

Current-main compatibility validation is recorded separately. The upstream rename to
CertifiedGaps is not backdated into any historical spectrum or validation metadata.
