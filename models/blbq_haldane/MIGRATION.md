# Core rename and historical evidence

This branch is based on current main, where the shared package is named `CertifiedGaps`
(UUID unchanged). Imports and the certificate Project.toml were migrated only inside this directory, which
was first developed as `certificates/blbq_haldane/` and later moved to `models/blbq_haldane/`
to follow the repository layout.

`migration_audit.jl` compares the old core at ef0a3e6 against current core sources after
normalizing the package name. It asserts equality for recurrence, ED, seed, RP, spin and
results modules, and checks that every noncore Manifest entry is unchanged. Its output is
`results/migration_audit.json`. Actual tests against the renamed package are recorded in
`results/validation_current_main.json`, separately from historical validation.

Original complete-spectrum/reference metadata are not rewritten to claim use of the
new name. The original validated source files were recovered and SHA-256 checked before
saving under `historical/thermal_source/`; the initial validation and analysis JSON are
preserved unchanged alongside them. Lossless gzip storage conversion has its own hash
record in `results/archive.json`.

One pilot's source file was optimized while its already-loaded Householder implementation
was running. Its raw end-of-run metadata is preserved under historical/. The result's
loaded-source hash is explicitly corrected to the archived Householder source, whose
hash also matches the independently completed L12 pilot. This correction and the observed
post-run disk hash are recorded in that result; no historical hash is silently replaced.
New long-running sector jobs capture their source hashes before starting computation.
