# Hexagonal AKLT: thermal-recursion pilot (no certificate)

`AKLTPilot` tests the 2+1-dimensional thermal recursion of
[the 2+1D note](../../outputs/thermal-mpo-gap-2plus1d.tex) on the periodic spin-3/2 AKLT model
of the honeycomb lattice, \(H=\sum_{\text{bonds}}P_{ij}(\text{total spin}=3)\), with this
projector normalisation throughout.

## Model

Each primitive cell contains `A(x,y)` and `B(x,y)`; the neighbours of `A(x,y)` are `B(x,y)`,
`B(x-1,y)` and `B(x,y-1)`, both coordinates periodic. The pilot uses 2×2, 4×2 and 4×4 cells (8,
16 and 32 spins). These small periodic quotients have short winding cycles, so they are
implementation checks and finite-size diagnostics, not bulk evidence. The three bond directions
are perfect matchings; a snake ordering is used only for the MPO, and every periodic bond is
retained.

The thermal code (`src/thermal.jl`, on the core's `purified_log_partition`) evolves the MPO
amplitude \(e^{-\beta H/2}\) from the normalised identity, applying gates only to the ket
indices. Each bond gate is the exact identity \(e^{-tP}=I+\operatorname{expm1}(-t)P\); the
noncommuting product still has second-order Trotter error and truncation error. Accumulated
Hilbert–Schmidt norms reconstruct the full, all-sector \(\log Z\); values with \(\log Z<0\) are
flagged, not clipped. No ground energy is fitted or subtracted (the exact AKLT ground energy is 0).

The reference diagonalisation (`src/exact.jl`) builds the projector from its polynomial in
\(\mathbf S_i\cdot\mathbf S_j\), independently of the Clebsch–Gordan form used by the gates, and
covers every total-\(S^z\) sector, counting \(S^z<0\) through spin flip. "Full" means no
Hilbert-space truncation; eigenvalues are still `Float64`.

For base geometry \((n,m)\) and inverse temperature \(\beta\) the diagnostics are
\(R=Z(n,m,\beta)-1\) (weight outside one zero-energy state), the doubling coefficient
\(C=[Z(2n,2m,2\beta)-1]/R^2\), the gluing defect \(G=\log Z(2n,2m,2\beta)-4\log Z(n,m,2\beta)\)
split as \(G=2G_x+G_y\), and \(\max(G,0)/R^2\). The relevant condition is \(CR<1\) with a uniform
bound on later scales; a small \(C\) is useless when \(R\) is large.

## Result

- **No gap certificate and no certified partition-function bound.** Every number this package
  produces is a `Float64` computation, stored with `certified=false`.
- **Negative result: the thermal-MPO purification does not resolve the smallest reference.** At
  \(\beta=2\) on eight spins the exact \(\log Z=5.31205522615\); with \(\delta\beta=1/8\), bond
  dimensions 16, 32, 64 and 128 miss it by about \(-1.773\), \(-1.424\), \(-0.587\) and
  \(-0.333\). Halving the time step does not remove the discrepancy, so it is compression bias,
  not Trotter error. Keeping moved sites between nonlocal gates changes the compression path but
  not the conclusion (\(\log Z=4.91941051889\) versus \(4.97859393485\) for the default at
  \(\chi=128\)). The scan temperatures are in any case far above the seed regime.
- **The gluing diagnostics are therefore unresolved.** Negative gluing defects \(G\) or small
  doubling coefficients \(C\) seen at one bond dimension must not be read as evidence for an
  all-scale inequality.
- **Eight-spin seed candidate (floating point).** The full 65,536-level spectrum of the 2×2-cell
  quotient has one numerically zero ground energy (\(1.35\times10^{-14}\)) and gap
  \(0.0920298633241\); its thermal excess is \(R(64)=0.00831113960\), and \(R=0.01\) is reached
  near \(\beta\approx61.9932\). This is an unverified dense eigensolve on a small quotient, not a
  certificate or bulk evidence.
- **Width-two projection-filter transfer (local Python, not reproducible from this repository).**
  The note's width-two analysis rests on two local Python checks that are not part of the
  repository (it has no Python): one finds that the one-step projection-filter spatial transfer
  (dimension 256) admits no positive metric for \(T\) or \(T^2\), so the self-adjoint transfer
  route fails for this filter; the other bounds the transfer spectrum with exact characteristic
  polynomials and gives \(d_{2L}\le d_L^2/4\) for every even longitudinal length \(L\ge4\) at
  width 2. Neither yields a bulk gap, and neither has been ported to Julia. The note's statement
  that these scripts live in `models/aklt/scripts/` is out of date.

The next requirement is a contraction method that reproduces the eight-spin reference and keeps
that accuracy as the width grows; until then this directory is an implementation check.

## Reproduce

`AKLTPilot` uses the shared core `CertifiedGaps` (`../../julia`, a path dependency) for spin
matrices, sector diagonalisation, purification and result output, with ITensors/ITensorMPS
pinned to the same versions as `models/spin2`. From the repository root:

```sh
julia --startup-file=no --project=models/aklt -e 'using Pkg; Pkg.instantiate()'
julia --startup-file=no --project=models/aklt models/aklt/test/runtests.jl   # 38 tests, ~3.5 min
julia --startup-file=no --project=models/aklt models/aklt/scripts/run.jl \
  models/aklt/configs/exact-small.toml models/aklt/results/exact-small.json
```

The other configurations are `thermal-pilot.toml` (2×2 base, \(\chi\in\{16,32\}\),
\(\delta\beta\in\{1/8,1/16\}\)), `thermal-reference-refinement.toml` (eight spins only,
\(\chi\in\{64,128\}\)) and `thermal-swap-refinement.toml` (the same without moving sites back).
Each output records source and environment hashes, parameters, versions, thread counts,
geometry, a `complete` flag and `certified=false`, and is rewritten after every finished sector
or temperature, so an interrupted run is visible; use a fresh filename to keep an earlier run.

## Layout

- `src/AKLTPilot.jl` (module), `src/model.jl` (periodic honeycomb quotient geometry),
  `src/thermal.jl` (thermal MPO purification), `src/exact.jl` (polynomial projector and
  reference diagonalisation).
- `scripts/run.jl`: runs one configuration from `configs/*.toml` and writes its JSON output.
- `test/runtests.jl`: the tests.

`results/` is gitignored and kept locally. The stored files were produced before the move to
`models/`, with the `spin2` environment, so their `environment/` hashes refer to
`spin2/Project.toml` and `spin2/Manifest.toml`; `pilot-summary.json` and the two projection
files there came from the local Python helpers mentioned above. Rerun the configurations to
regenerate the Julia outputs.
