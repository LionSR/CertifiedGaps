# Draft PR publication and audit status

The original completed computation is preserved at commit
`6c375b5981f1f7c7813bce53ca01a1ce1c04c9f2`. Its core was named ThermalMPOGap.
Upstream main subsequently renamed the same package UUID to CertifiedGaps.
GitHub rejected the initial push under its private-email protection. Unpublished
commit identities were changed to the authenticated account’s GitHub noreply
address. The archived computation commit above has exactly the same tree
(`7ec89f1635d7b9a493159bd8433b7459f93bb3bf`) as the original local commit
`107116f7b5cc79db5023553f14dec2b98245d375`; no historical file or result hash
was changed. Only the current checker’s archive locator was updated.

This branch incorporates upstream main without proposing core changes, and its
loader aliases CertifiedGaps under the historical name. The original JSON files
and their source hashes are preserved unchanged. The fast checker verifies those
hashes against the archived original commit, then recomputes the RP identities,
seed conversions and recurrence using the current core. Thus a full Git history
(or access to that commit) is required for the historical-source audit.

The independent [audit](AUDIT.md) found no substantive mathematical blocker.
It supplies explicit cell Gram, odd-seed recurrence/interpolation and uniqueness
arguments, and reviews symmetry counts and enclosure arithmetic. The fresh
twelve-spin spectral replay completed successfully in 949.35 seconds (15.82 minutes).
The replay recomputes the model matrices and dense eigensystems; it does not use
the stored eigenvalue radii or partition enclosures. It preserves original outputs:

```sh
OPENBLAS_NUM_THREADS=1 \
SPIN1_DIMERIZED_REPLAY_PATH=certificates/spin1_dimerized/results/replay_certificate_d7.json \
julia --project=julia certificates/spin1_dimerized/certify.jl 7
```

The stored output reports a lower bound of 0.2187941698399417 at delta=7/10,
n=3, beta=17/4 for even cell lengths L>=3, i.e. physical lengths 8,12,16,... .
The PR remains a draft for human review, with the independent audit and replay
now recorded. The odd-cell extension
is a mathematical argument in this directory, not a separately machine-checked
Lean theorem. No Lean formalization of these numerical inputs is claimed.

Audit limitations: the fast checker accepts stored log-partition enclosures as
inputs, verifies their rational seed conversion, and reruns the recurrence; it
does not repeat the large spectral enclosure computation. Numerical reference
gaps, finite-grid failures and finite-size extrapolations are not certificates.
The explicit rounding argument assumes ordinary IEEE binary64 arithmetic and
classical BLAS dot products; these assumptions are part of the reviewer audit.
All dense spectra are retained (about 73 MB of JSON) for inspection.

Only this certificate directory is proposed. No core changes, CI additions,
merges into the target branch, or ready-for-review promotion are authorized here.

## Current-core compatibility validation

Against upstream `22b2759e0790197b83e2571b465aa8dc69ebc0be`, the main checker
passed 115 assertions and the supplemental verified-block checker passed 16.
[The separate validation record](results/current_main_validation.json) preserves
current source hashes and log fingerprints. The core recurrence implementation
is unchanged from the archived computation; relevant module-header differences
are documentation names. Certificate arithmetic is unchanged. The loader, local
package dependency name, historical-source lookup and optional replay destination
are the compatibility changes. The user-requested `certificates/spin1_dimerized/`
location is retained despite the new general `models/` convention.

## Completed full replay

The separate [replay result](results/replay_certificate_d7.json) records a fresh
dense eigensystem and enclosure calculation for all 7 small-ring and 79
large-ring blocks, with weighted dimensions 729 and 531441.
`replay_check.jl` passed all 40 assertions in its 2.0-second test set: every
mathematical output (including every block error bound and log-partition endpoint)
is identical to the original, the current source hashes match, and the rational
seeds reproduce the bound. Only provenance and runtimes differ.
[The comparison record](results/replay_comparison.json) stores both result hashes
and runtimes. Original result files remain unchanged.

```sh
OPENBLAS_NUM_THREADS=1 julia --project=julia certificates/spin1_dimerized/replay_check.jl
```

There are no remaining execution or audit blockers. This is a computational
certificate supported by an explicit mathematical proof and independent review,
subject to the stated arithmetic assumptions; it is not a Lean-checked model
instantiation. No changes have been merged into the target branch.
