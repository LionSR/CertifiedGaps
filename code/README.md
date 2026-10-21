# Code

| Directory | Purpose |
|---|---|
| [`julia/`](julia/) | Shared `CertifiedGaps` Julia package and tests |
| [`models/`](models/) | Model-specific generators, checkers, tests, and method documentation |
| [`lean/`](lean/) | Criterion proofs and kernel-checked seed arithmetic |
| [`leaderboard/`](leaderboard/) | Lean gap records and the theorem audit |
| [`tools/`](tools/) | Upstream evidence fetcher and leaderboard validation |

Run commands from the repository root; see [environment setup](../environment/README.md)
and each model's README.

The Julia and Lean source trees retain their relative positions. This preserves model manifests,
cross-model imports, and the source paths embedded in existing certificates. Historical source
comments and the Lean README use paths relative to this `code/` workspace.

Each `models/<slug>/results` is a relative symlink to `../data/models/<slug>/` from this
directory. The canonical files are in `data/`; the links let hash-pinned generators and checkers
continue using their original model-local paths. Likewise, `references` links to
`../data/references/`. Preserve these symlinks when copying or checking out the repository.
