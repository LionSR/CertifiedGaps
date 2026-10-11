# Validation against CertifiedGaps

The draft branch is based on main commit
`22b2759e0790197b83e2571b465aa8dc69ebc0be`. Its `AGENTS.md`, `julia/README.md`
and package files specify `CertifiedGaps` (same UUID and version 0.2.0).
The user explicitly requires this certificate to remain under `certificates/`,
overriding main's new default `models/` destination. No shared files were edited.

## Preserved evidence

The old validated snapshot remains commit
`1c4c07711072b8bfe1a60d188a498a491b8dc42b` in the draft PR's history.
Every pre-migration result JSON and log, including the five certificates,
all spectral replay attempts, and `audit_manifest.json`, is retained byte-for-byte.
The old manifest describes the old snapshot; it is not a checksum claim about
the migrated active scripts. `historical_sources/` preserves the 47 executable
source witnesses named by each original certificate, including the old core.
These files are evidence, not a second active package or a core modification.

`historical_provenance.jl` verifies all 235 historical hash entries and asserts
that each corresponding current source is exactly the old bytes with
`ThermalMPOGap` replaced by `CertifiedGaps`. This includes the model code,
enclosure arithmetic, recurrence, original tests and shared source dependencies.
The three recurrence implementation files are byte-identical. Therefore no
mathematical algorithm changed in this migration. Upstream Manifest changes
rename the package and its extension; dependency versions are unchanged.

## Separate current validation

Run from the repository root with Julia 1.12.4 and `OPENBLAS_NUM_THREADS=1`:

```sh
julia --project=julia certificates/heisenberg_staggered_field/environment.jl
julia --project=julia certificates/heisenberg_staggered_field/test/runtests.jl
julia --project=julia certificates/heisenberg_staggered_field/audit_bounds.jl
julia --project=julia certificates/heisenberg_staggered_field/spectral_replay.jl
```

New records use the `migration_` prefix. The environment record contains the
base revision, historical revision, current source hashes and package versions.
The spectral replay reconstructs every block and recomputes the full spectral
certificate for all five fields without a scan or resumed historical witnesses.
The independent reviewer checks this migration separately. It remains the same
reviewed spectral implementation, not a second implementation.

The historical `spectral_replay_attempt1.jl` and `spectral_replay_attempt2.jl`
are left unchanged as failed-attempt evidence. Reproduce those only from the
old snapshot. The active replay's `--resume-attempt2` option checks archived
sources before using those old witnesses; the migration validation does not use
that option. `seal.jl` creates a separate `migration_manifest.json`.

## Observed outcome

All 71 fast certificate tests and 17 independent accumulation/normalization
tests passed under CertifiedGaps. The full five-field, 270-block fresh replay
completed successfully in 189.418 seconds, without resuming any old witness.
All five fresh raw bounds exactly matched the preceding completed replay and
support the unchanged six-decimal published bounds. The independent reviewer
inspected the migration source changes, provenance checker and completed replay
and accepted rename-only mathematical equivalence and all five bounds; the
reviewer did not execute a further migration replay.

The environment and migration checker verify 235 historical source entries.
The migration checker repropagates the fresh enclosures and tests the five
published claims; its saved JSON records the replay hash. All 106 pre-migration
result files were byte-compared against the old Git snapshot; 105 are covered
by the original manifest (the remaining file is that manifest itself).
New results/logs and the current deliverable manifest are separate. Run the
additional acceptance check with:

```sh
julia --project=julia certificates/heisenberg_staggered_field/migration_check.jl
```
