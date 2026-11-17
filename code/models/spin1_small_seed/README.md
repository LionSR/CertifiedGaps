# Spin-1 Heisenberg ring: a small seed from rings of at most 12 sites

This model certifies the leaderboard challenge *a small seed for spin-1*: the seed calculation
uses only twisted thermal traces on rings of at most 12 sites (even lengths, the two real twists
$`I`$ and $`P`$, absolute error $`10^{-3}`$), and certifies $`\Delta_L > 0.039036761`$ for every
even $`L \ge 60`$. Every scalar step is replayed in exact rational arithmetic from hash-audited
frozen spin-1 inputs, so the whole finite seed is small enough to be kernel-checked in Lean 4
(the Lean side remains open; see [`code/models/lean_seeds`](../lean_seeds/README.md)).

| parameter | seed (n₀, β₀) | certified gap ≥ | lengths covered |
|---|---|---:|---|
| periodic ring | $`(60, 49/2)`$ | $`0.039036761`$ | even $`L \ge 60`$ |

## Model

The periodic spin-1 antiferromagnetic Heisenberg ring of $`L`$ sites (lengths count sites),

```math
H_L = \sum_{j=1}^{L} \mathbf{S}_j \cdot \mathbf{S}_{j+1},
\qquad \mathbf{S}_{L+1} \equiv \mathbf{S}_1.
```

with unit coupling and no on-site term. The finite inputs are twisted thermal traces on short
rings, $`Z_n(\beta, g) = \mathrm{Tr}(g\, e^{-\beta(H_n + n\kappa I)})`$ with the rational
per-site shift $`\kappa = 700741/500000`$ and a global spin rotation $`g \in \{I, P\}`$ ($`P`$
the $`\pi`$ rotation about the $`y`$ axis) at the closing bond. The extensive shift cancels
identically in both purities, $`P_{\mathrm{s}}(n, \beta) = Z_{2n}(\beta)/Z_n(\beta)^2`$ and
$`P_\tau(2n, \beta) = Z_{2n}(2\beta)/Z_{2n}(\beta)^2`$.

## Reflection positivity

With the integer spin-1 matrices of `Spins.integer_spin_one` ($`\mathtt{xp} = \sqrt{2} S^x`$,
$`\mathtt{ap} = \sqrt{2} i S^y`$, $`\mathtt{sz} = S^z`$ in the $`S^z`$ eigenbasis), the bond is
$`h = \mathtt{sz} \otimes \mathtt{sz} + \tfrac{1}{2}(\mathtt{xp} \otimes \mathtt{xp} -
\mathtt{ap} \otimes \mathtt{ap})`$. The real symmetric involution $`R = \mathrm{antidiag}(1,
-1, 1)`$ reverses $`S^z`$ and $`S^x`$ and keeps $`S^y`$, giving the exact factorization

```math
-h = \mathtt{sz} \otimes R\,\mathtt{sz}\,R + \tfrac{1}{2}\,\mathtt{xp} \otimes R\,\mathtt{xp}\,R
     + \tfrac{1}{2}\,\mathtt{ap} \otimes R\,\mathtt{ap}\,R,
```

with on-site term $`o = 0`$. `rp_certificate` from the core verifies this identity exactly over
`Complex{Rational{BigInt}}` and throws otherwise; `check.jl` re-proves it on every run. The
Cartesian form $`-h = \sum_\alpha G_\alpha \otimes G_\alpha`$ is formalized in
`code/lean/CertifiedGaps/Heisenberg/Main.lean`.

## Seeds and certified bounds

The ten trace rows used (lengths $`n \in \{4, 6, 8, 10, 12\}`$, twists $`I, P`$,
$`\beta = 49/4`$) are read from the frozen spin-1 certificate `data/models/spin1/source1.json`,
whose SHA-256 is pinned in its `verification.json`. Each row carries a rational interval
$`[\mathrm{lower}, \mathrm{upper}]`$ around its center with half-width $`\le 3 \times
10^{-5}`$ (its recorded `allowed_error`), well inside the challenge tolerance $`10^{-3}`$, and
`matches_upstream_exactly = true`. The stored length-120 integer MPS trial
(`data/models/spin1/trials.json`, likewise hash-pinned) satisfies
$`500000\,U + 700741\,V < 0`$ with $`V > 0`$, i.e. energy per bond below $`-\kappa`$.

From the two-twist sector inversion $`J_n = (Z_n + 3P_n)/4`$, $`N_n = (Z_n - P_n)/4`$ with
per-sector errors $`\varepsilon \le 10^{-3}`$, degree-4 squared filters $`F(x) = x^4 q(x)^2`$
with $`q_J(x) = x^4 - 0.32529 x^2 + 0.008718`$ and $`q_N(x) = x^4 - 0.386904 x^2 + 0.011397`$
(exact rationals; every exterior ray is re-proved exactly), and the temperature lift
$`49/4 \to 49/2`$ at moment $`k = 30`$ with the length-120 trial, the exact seed at
$`(n_0, \beta_0) = (60, 49/2)`$ has defects

| quantity | exact value (rounded display) |
|---|---|
| $`\varepsilon_{\mathrm{s}} = p \ge 1 - P_{\mathrm{s}}(60, 49/2)`$ | $`0.047914767501\ldots`$ |
| $`\varepsilon_\tau = q \ge 1 - P_\tau(120, 49/2)`$ | $`0.112578547343\ldots`$ |

both below the diagonal threshold $`\varepsilon_* \approx 0.12731`$. The 9-step rational
doubling chain on the $`10^{-512}`$ rounding grid with rate $`r = 130115651/50000000`$ at
$`\beta_0 = 49/2`$ then certifies

```math
\Delta_L > \frac{\log r}{\beta_0} > \frac{39036761}{10^9} = 0.039036761
\qquad \text{for every even } L \ge 60,
```

and dominates the published target $`\log(20)/784`$ (checked as the exact integer inequality
$`r^{32} > 20`$). The full record is `results/certificate.json` (exact rational strings, source
hashes, and run metadata).

## Reproduce

From the repository root with Julia 1.12 and `OPENBLAS_NUM_THREADS=1` (times measured on an
8-core x86-64 workstation, 2026-10-30):

```sh
julia --project=code/julia code/models/spin1_small_seed/certify.jl            # 19 s
julia --project=code/julia code/models/spin1_small_seed/check.jl                # 17 s
julia --project=code/julia code/models/spin1_small_seed/check.jl --full         # 18 s
julia --project=code/julia code/models/spin1_small_seed/test/runtests.jl        # 23 s, 26 tests
```

`check.jl` re-audits the frozen inputs, re-proves the RP identity, replays the exact chain,
and requires every stored field of `results/certificate.json` to match the recomputation;
`--full` additionally audits all 18 source-1 trace rows. Any mismatch throws.

## What is rigorous

- **Trusted base: exact rational arithmetic throughout.** The twisted-trace intervals and the
  MPS trial energy are inherited from the frozen, hash-audited spin-1 certificate
  (`data/models/spin1/source1.json`, `trials.json`; component SHA-256 pinned in its
  `verification.json` and re-verified here), where they were produced by exact integer
  arithmetic with explicit a priori rounding bounds and checked against the upstream preprint's
  printed integers. A from-scratch recomputation of the traces is
  `code/models/spin1/check.jl --full` (about 45 min); it is the same frozen trusted base this
  model audits.
- **Exact in this model:** the reflection-positivity identity (`rp_certificate`), the sector
  inversion, the filter exterior rays, the capped masses, the temperature lift, the seed
  defects, and the doubling chain — all `Rational{BigInt}`, no floating-point comparison.
- **Floating-point only:** nothing in the certificate. The filter coefficients were proposed by
  a `Float64` search in `code/models/lean_seeds` (design C) and are re-proved exactly here.
- **Platform caveats:** none; no BLAS or LAPACK is used, so results are bit-identical across
  platforms.

## Layout

Stored results live in [`data/models/spin1_small_seed/`](../../../data/models/spin1_small_seed/);
this directory's `results/` is a relative symlink to that canonical location.

| path | role |
|---|---|
| `certify.jl` | regenerates `results/certificate.json` from the frozen inputs |
| `check.jl` | verifies the certificate (`--full` also audits all 18 source-1 trace rows) |
| `src/Spin1SmallSeed.jl` | the exact certificate module (RP, input audit, chain, writer/checker) |
| `test/runtests.jl` | 26 unit and regression tests |
| `IDEAS.md` | route memo for the challenge |
| `results/certificate.json` | the certificate: seed, defects, chain, gap bound, hashes, metadata |

No hash-pinned files beyond `results/certificate.json`; no own `Project.toml` (the model runs
in the core environment `--project=code/julia`).

## History

- 2026-10-30: created for the leaderboard challenge `spin1-small-seed`, replaying design C
  (`real_even`) of `code/models/lean_seeds` as a standalone certificate with its own checker
  and leaderboard record.
