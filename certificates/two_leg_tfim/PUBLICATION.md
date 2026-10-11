# Draft PR integration status — 2026-10-11

The scientific result and independent replay are complete at the audited revisions:
certificate `d546444493b82f8bd9990d55dec585c2bb3ef30e`, replay audit
`7f03adf4217e4947c6a73b7b3ed07dbede079629`. Their baseline is
`ef0a3e6`, with the Julia package named `ThermalMPOGap`.

At publication, upstream `main` is `22b2759e0790197b83e2571b465aa8dc69ebc0be`.
That revision renames the core to `CertifiedGaps` and changes core/project/source
fingerprints. The draft deliberately retains the verified branch and historical
artifacts. Its scripts and strict source-hash checker have **not** been migrated or
validated against the renamed core. Merging as-is would leave imports incompatible
with current main. This is an integration blocker, not a failure of the verified
bound on the audited branch.

Before making the PR ready: adapt imports inside this certificate directory to the
renamed package, distinguish historical provenance from current source checks,
validate the migration, and record any new replay separately. Do not rewrite old
certificate/audit JSONs or change shared core files merely to hide hash differences.

For exact reproduction now, check out this PR's head or the audited commits and
follow `README.md` and `AUDIT-2026-10-11.md`. No active computation remains. The PR
changes only `certificates/two_leg_tfim/`; it adds no CI and makes no core edits.

## Privacy-preserving publication commit mapping

GitHub rejected the initial push under its private-email protection. The publication
branch replaces only author/committer email metadata with the authenticated account's
GitHub noreply address. Original commits remain on local branch
`audit/two-leg-tfim-original`; historical JSONs and audit notes are unchanged.
The following published commits have **identical Git trees** to their originals:

| Historical commit | Published equivalent |
|---|---|
| `d546444493b82f8bd9990d55dec585c2bb3ef30e` | `bf7efa7711b317aad3ec04a40c728c075083d1a0` |
| `7f03adf4217e4947c6a73b7b3ed07dbede079629` | `8a1aad63e1a285392c97bf5cf17e8da03f0b8314` |
| `73a22841f28f4952ab0589d988665f0bcd9cecc1` | `f776c3d13f75f6c4cc0be69985f8d4412f480b01` |

When reproducing from GitHub, check out the PR head or these published equivalents.
No data, mathematical implementation, or stored fingerprint changed in this rewrite.
