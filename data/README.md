# Data

Stored results have one canonical location, `data/models/<slug>/`. Each corresponding
`code/models/<slug>/results` is a relative symlink so existing generators and checkers can
read the same files without changing their source hashes. The reorganization does not alter
certificate contents or recompute numerical results.

| Dataset | Origin and use |
|---|---|
| `models/spin1/` | Exact trace, trial-state, and scalar data reproduced from the pinned OpenAI evidence; checked by `code/models/spin1/check.jl` |
| `models/heisenberg_staggered_field/`, `models/dimerized_heisenberg/`, `models/spin1_dimerized/`, `models/spin1_largeD/`, `models/heisenberg_ladder/`, `models/two_leg_tfim/`, `models/tfim_longitudinal/` | Model generators produce seed enclosures and certificates; each model's `check.jl` replays its certificate arithmetic |
| `models/lean_seeds/` | Witness studies, exact arithmetic outputs, and cost diagnostics supporting the Lean seed generators |
| `models/spin2/` | Stored trial certificates and local numerical studies; no uniform spin-2 gap certificate |
| `models/blbq_haldane/` | Negative-result scan summaries and diagnostic data; no gap certificate |
| `models/aklt/` | Local, regenerable pilot outputs; no gap certificate and no shipped numerical dataset |
| `references/` | External evidence and paper sources; see below |

Each [model README](../code/models/) identifies its data, trusted arithmetic, generation commands,
checker commands, and measured runtimes. `check.jl --full` replays spectra where supported.
Diagnostic scans are not certificates. Regenerable outputs are ignored by `data/.gitignore`.

Spin-1 experiment datasets, including `temperature-lift.json`, live under
`data/models/spin1/experiments/`. Relative file symlinks beside the scripts preserve the
paths used by their existing consumers.

## External evidence

The upstream source is
[`openai/math` at `adc7f1241b42e322a6451854ab7e4b4c146bf78a`](https://github.com/openai/math/tree/adc7f1241b42e322a6451854ab7e4b4c146bf78a).
From the repository root:

```sh
# Required by the spin-1 default checker; hashes are verified by the fetcher.
julia code/tools/fetch_upstream.jl

# Additional sources for the Lean seed studies.
julia code/tools/fetch_upstream.jl --boundary
julia code/tools/fetch_upstream.jl --paper
```

The first command downloads periodic-gap evidence to `data/references/periodic-gap-source/`.
`--boundary` downloads trial data to `data/references/boundary-field-gap-source/`.
`--paper` writes the selected upstream LaTeX sources to `code/models/lean_seeds/upstream/`.
All downloads stay ignored. Other models' default checks do not require these downloads.
