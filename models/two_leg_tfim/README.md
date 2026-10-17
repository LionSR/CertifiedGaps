# Two-leg transverse-field Ising ladder

This module certifies a uniform spectral gap for the two-leg transverse-field Ising ladder at
transverse field $g = 3$: **$\Delta_L \ge 0.636870831$ and a unique ground state for every even
$L \ge 4$ rungs** ($N = 2L$ spins, $N$ divisible by 4, $N \ge 8$). The certificate is produced
end to end by the shared core (`CertifiedGaps.GapCertificates.certify_gap`) from the ladder
written as a reflection-positive `RingModel` of rungs. The method, scan tables and finite-size
diagnostics are in [`METHOD.md`](METHOD.md); the route analysis is in [`IDEAS.md`](IDEAS.md).

| $g$ | seed $(n_0 \text{ rungs}, \beta_0)$ | certified gap $\ge$ | lengths covered |
|---|---|---:|---|
| 3 | $(4, 3/2)$ | 0.636870831 | every even $L \ge 4$ rungs |

## Model

Each unit cell is a two-spin rung ($d = 4$). In Pauli normalization on a periodic ring of $L$
rungs,

```math
H_L = -\sum_{j=1}^L (Z_{j,a} Z_{j+1,a} + Z_{j,b} Z_{j+1,b} + Z_{j,a} Z_{j,b}) - g \sum_{j=1}^L (X_{j,a} + X_{j,b}).
```

Lengths count rungs, and the $L = 2$ ring retains both periodic leg bonds, following the core
convention. The certificate works in the **rotated frame** obtained by a Hadamard gate on every
spin, i.e. the exact rational orthogonal involution $W = (H \otimes H)/2$ with
$`H = \begin{pmatrix}1&1\\1&-1\end{pmatrix}`$ on every rung (`Ladder.rung_hadamard`). It maps
$X \leftrightarrow Z$, so the rotated ladder (`Ladder.ladder(g)`) has bond
$h = -(X_a X_{a'} + X_b X_{b'})$ and on-site term $o = -X_a X_b - g(Z_a + Z_b)$. There the
$`\mathbb{Z}_2`$ charge (number of down spins mod 2, `charge = [0, 1, 1, 0]`) is conserved and
the leg exchange $`a \leftrightarrow b`$ is an exact on-site $`\mathbb{Z}_2`$ symmetry
(`symmetries = [LEG_SWAP]`). The tests check $`W o_{\text{phys}} W = o_{\text{rot}}`$ and
$(W \otimes W) h_{\text{phys}} (W \otimes W) = h_{\text{rot}}$ exactly. $W^{\otimes L}$ is
orthogonal, so both frames have the same spectrum and purities.

## Reflection positivity

In the rotated frame, the real bond factors $C_1 = X \otimes I$ and $C_2 = I \otimes X$ give

```math
-h = C_1 \otimes C_1 + C_2 \otimes C_2,
```

and the on-site term $o$ is real symmetric; a single rung needs no internal reflection ($F = I$).
`Ladder.ladder_rp` passes these factors to the core `rp_certificate`, whose `verify_rp` checks
the identity over `Rational{BigInt}`; the certificate stores them and `check.jl` re-verifies
them. In the physical frame the same identity holds with $C_1 = Z \otimes I$, $C_2 = I \otimes Z$
(the form used in `lean/CertifiedGaps/Models/TwoLegIsing.lean`), which the tests also check.

## Seeds and certified bounds

| $g$ | seed $(n_0, \beta_0)$ | $\varepsilon_{\mathrm{s}} \ge 1 - P_{\mathrm{s}}(4, 3/2)$ | $\varepsilon_\tau \ge 1 - P_\tau(8, 3/2)$ | certified gap $\ge$ (interval) | exact chain bound | 8-rung ED gap |
|---|---|---:|---:|---:|---:|---:|
| 3 | $(4, 3/2)$ | $0.11238966699069$ | $0.051536220063797$ | 0.636870831 | $> 0.636342881$ | 2.6618 (numerical) |

The seeds are the exact rationals `seed.eps_s = 11238966699069//100000000000000` and
`seed.eps_t = 51536220063797//1000000000000000` of
[`results/certificate_g3_1.json`](results/certificate_g3_1.json). The headline is the
outward-rounded interval bound `certified_gap_bound(εₛ, εₜ, 3//2).bound` $= 0.6368708319575064$
(`interval_bound_decimal`, rounded down above), about $`24\%`$ of the 8-rung numerical gap; the
certificate's claim `gap_lower_bound` is the exact rational doubling chain,
$\Delta_L > 636342881/10^9$. Both are rigorous. The defects satisfy the Lean theorem's seed
hypotheses ($\varepsilon_{\mathrm{s}} \le 0.1125$, $\varepsilon_\tau \le 0.05162$, checked
exactly in the tests). On the scanned grid $`g \in \{2, 5/2, 3, 4\}`$, the fields $g = 2$ and
$5/2$ have no admissible seed with $n_0 \le 4$ rungs, and $g = 4$ (numerical bound $2.142$) was
not certified because $g = 3$ already is.

## Reproduce

From the repository root, with Julia 1.12 and `export OPENBLAS_NUM_THREADS=1`, in the core
environment:

```sh
julia --project=julia -e 'using Pkg; Pkg.instantiate()'
julia --project=julia models/two_leg_tfim/certify.jl          # results/certificate_g3_1.json
julia --project=julia models/two_leg_tfim/check.jl            # fast checker
julia --project=julia models/two_leg_tfim/check.jl --full     # + freshly verified spectra
julia --project=julia models/two_leg_tfim/test/runtests.jl    # same as the fast checker
# Numerical seed scan for one field, then the extrapolation and the METHOD.md tables:
julia --project=julia models/two_leg_tfim/scripts/scan.jl 3//1
julia --project=julia models/two_leg_tfim/scripts/extrapolate.jl
julia --project=julia models/two_leg_tfim/scripts/summarize.jl
```

`certify.jl [g]` certifies every point of `Ladder.POINTS` (or the one at field `g`). The scan
writes full spectra to `results/spectra_g<g>.json` (gitignored) for `extrapolate.jl` and
`summarize.jl`. Measured run times (Julia 1.12.4, one BLAS thread, Apple M-series laptop shared
with other jobs):
- `certify.jl`: 532 s and 980 s in two identical runs at load averages of about 50–80 (32
  blocks; the 8-rung blocks have dimension at most 2102, most of them complex);
- `check.jl` / `test/runtests.jl`: 37 tests, 38–45 s including compilation;
- `check.jl --full`: 319 s (284 s to recompute both verified spectra);
- `scripts/scan.jl`: about 5 min per field (unverified spectra of the core blocks).

## What is rigorous

- **Rigorous:** the exact RP identity and frame change (`Rational{BigInt}`); the complete
  spectra of the rings of 4 and 8 rungs, split exactly by the core `VerifiedSpectra.ring_blocks`
  into $`\mathbb{Z}_2`$ charge, lattice momentum (±k merged with multiplicity 2) and leg-exchange
  parity, with block entries enclosed in interval arithmetic; each block enclosed by
  `CertifiedGaps.VerifiedSpectra.verified_eigen` (residual and Gram matrices in `Float64` with
  a priori $`\gamma_{N+8}`$ (real blocks) or $\gamma_{2N+8}$ (complex blocks) rounding bounds and
  Frobenius norms, then Ostrowski/polar-factor and Weyl bounds, plus the exact construction
  error); $\log Z$ by `log_partition_enclosure` in interval arithmetic; seed defects rounded up
  to a $10^{-15}$ grid; the exact rational doubling chain and the interval
  `certified_gap_bound`. Trusted base: IEEE binary64 round-to-nearest with gradual underflow,
  standard (non-fast-math) BLAS and LAPACK producing finite output, `IntervalArithmetic`, and
  the repository's two-purity theorem.
- **Fast vs. full verification:** `check.jl` re-verifies the RP identity, the stored model
  against `Ladder.ladder(g)`, the exact chain from the stored seeds, the interval bound, the
  Lean seed hypotheses and the model-local source hashes; `--full` also recomputes and
  re-verifies every block spectrum and requires the fresh defects to be at most the stored ones.
- **Numerical only (not used in the proof):** the seed search ($`g \in \{2, 5/2, 3, 4\}`$,
  $`n_0 \in \{2, 3, 4\}`$, $\beta = 1/8 : 1/8 : 32$), the negative scan results at $g = 2$ and
  $5/2$, the extrapolated seed sizes, finite-ring ED reference gaps, and the literature critical
  field $g_c \approx 1.83214(5)$.

## Layout

- `certify.jl`, `check.jl`: certificate generator (core `certify_gap`) and checker (core
  `check_certificate`).
- `src/Ladder.jl`: the model only: physical and rotated terms, the frame change, the core
  `RingModel` with its charge and leg-exchange symmetry, the RP factors, the certified points,
  the hashed source set, and unverified energies for scans.
- `scripts/`: `scan.jl` (numerical seed scan), `extrapolate.jl`, `summarize.jl` (rewrites the
  `## Results` section of `METHOD.md`).
- `test/runtests.jl`: exact frame and RP tests, core blocks against an independent full
  Hamiltonian, and the stored-certificate checks.
- `results/certificate_g3_1.json`: the core-format certificate (model, RP data, seed, log Z
  enclosures, spectrum summaries, exact chain, and `metadata` with the model-local
  `source_sha256` of `src/`, `certify.jl`, `check.jl` and the core `core_sha256`);
  `scan_g*.json` (numerical scans) and `extrapolation.json`.

## History

The original certificate (`results/certificate_g3_1.json.gz`, generated when the directory was
`certificates/two_leg_tfim`) used local symmetry blocks (parity, leg exchange, reflection, half
translation), a local polar/Weyl enclosure and MPFR log Z. It was retired together with its
replay scripts, `environment.json` and the model's own Julia environment. The model now runs
entirely on the shared core (`certify_gap` with the leg-exchange symmetry) and was regenerated
with the same seed: $\varepsilon_{\mathrm{s}}$ $0.1124084 \to 0.1123897$,
$\varepsilon_\tau$ $0.0516152 \to 0.0515362$, interval bound $0.636387 \to 0.636871$ (exact
chain bound $0.636343$).
