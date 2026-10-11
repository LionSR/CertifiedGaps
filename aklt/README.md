# Hexagonal AKLT thermal-recursion pilot

This pilot tests the thermal recursion in
[the 2+1D note](../outputs/thermal-mpo-gap-2plus1d.tex) for the periodic
spin-3/2 AKLT Hamiltonian `H = sum_bonds P(total spin = 3)`.
All energies use this projector normalization. Results are numerical,
not certified gap or partition-function bounds.

## Geometry and implementation

Each primitive cell contains `A(x,y)` and `B(x,y)`. The three neighbors
of `A(x,y)` are `B(x,y)`, `B(x-1,y)`, and `B(x,y-1)`, with both coordinates
periodic. The pilot uses 2x2, 4x2, and 4x4 cells: 8, 16, and 32 spins.
These are small periodic quotients; short winding cycles need not have
the same geometry as local hexagons in the infinite lattice. They are
implementation checks and finite-size diagnostics, not bulk evidence.

The three bond directions are perfect matchings. A snake ordering is
used only for the MPO; all physical periodic bonds are retained.

The thermal code evolves an MPO amplitude `exp(-beta H/2)` from the
normalized identity, applying gates only to its ket indices. Each bond
gate is evaluated by the exact projector identity
`exp(-t P) = I + expm1(-t) P`. The complete noncommuting evolution still
has second-order Trotter error and tensor truncation error. Accumulated
Hilbert-Schmidt norms reconstruct the full, all-sector `log Z`.
Values with `log Z < 0` are retained and flagged; they are not clipped.

The independent diagonalization assembles the bond projector from its
polynomial in `S_i . S_j`, rather than the Clebsch-Gordan formula used
by the tensor gates. It includes every total-Sz sector, using spin-flip
symmetry to count multiplicities. "Full diagonalization" means no
Hilbert-space truncation; eigenvalues remain Float64 approximations.

## Run

`aklt/` is the Julia package `AKLTPilot` with its own environment. It uses
the shared core `ThermalMPOGap` (`../julia`, a path dependency) for spin
matrices, sector diagonalization, thermal purification and result output,
with ITensors/ITensorMPS pinned to the same versions as `spin2/`. From the
repository root:

```sh
julia --startup-file=no --project=aklt -e 'using Pkg; Pkg.instantiate()'
julia --startup-file=no --project=aklt aklt/test/runtests.jl
julia --startup-file=no --project=aklt aklt/scripts/run.jl \
  aklt/configs/exact-small.toml aklt/results/exact-small.json
julia --startup-file=no --project=aklt aklt/scripts/run.jl \
  aklt/configs/thermal-pilot.toml aklt/results/thermal-pilot.json
```

Format only files you change, with `margin = 100` (JuliaFormatter is in the
`spin2` environment), e.g.
`julia --project=spin2 -e 'using JuliaFormatter; format_file("aklt/src/exact.jl"; margin = 100)'`.

The stored results were produced before the migration, with
`--project=spin2` and the pre-migration sources; their `environment/` hashes
refer to `spin2/Project.toml` and `spin2/Manifest.toml`.

Each output records source and environment hashes, parameters, versions,
thread counts, geometry, and `certified=false`. Results are saved after
each completed magnetization sector or temperature target. `complete`
distinguishes finished from interrupted runs. Use a fresh output filename
to preserve an earlier run.

## Quantities to inspect

For base geometry `(n,m)` and inverse temperature `beta`:

- `R = Z(n,m,beta)-1` measures total weight outside one zero-energy state.
- `C = [Z(2n,2m,2beta)-1]/R^2` is the actual finite-step doubling coefficient.
- `G = log Z(2n,2m,2beta)-4 log Z(n,m,2beta)` is the aggregate gluing defect.
- `G_x` and `G_y` resolve the two spatial gluing steps; `G = 2 G_x + G_y`.
- `eta_over_R_squared = max(G,0)/R^2` tests the scale of an approximate-gluing error.

The initial grid crosses two time steps with two bond dimensions, so
these two errors can be inspected independently. It is a short-temperature
implementation pilot. A small `C` at high temperature is not useful if
`R` is enormous: the relevant condition is `C R < 1`, and a proof would
need a uniform bound on subsequent scales. Negative gluing defects at
one bond dimension are not rigorous inequalities.

Full-trace normalization and preservation of small excited weights are
essential. Convergence of the ground-state energy alone does not validate
these measurements. The exact AKLT ground energy is zero; the simulation
does not fit or subtract a numerical ground-state energy.

## First results

The full eight-spin spectrum has 65,536 eigenvalues, one numerically
zero ground energy (`1.35e-14`), and gap `0.0920298633241`. Its thermal
excess is `R(64) = 0.00831113960`; the numerical `R = 0.01` threshold is
approximately beta `61.9932`. This is a finite-seed candidate, not a
certificate or a bulk-gap result.

The first thermal scans do **not** resolve the gluing criterion. At
beta 2 the eight-spin reference is `log Z = 5.31205522615`. With time
step 0.125, bond dimensions 16, 32, 64, and 128 give errors approximately
`-1.773`, `-1.424`, `-0.587`, and `-0.333`, respectively. Halving the
time step does not remove the discrepancy. The scan temperatures are
also far above the low-temperature seed regime.

A separate comparison keeps moved sites between nonlocal gates and
restores their order at the end of each step. This changes the
compression path while preserving untruncated evolution. It does not
fix the reference error: at beta 2 and bond dimension 128 it gives
`log Z = 4.91941051889`, versus `4.97859393485` for the default policy.

The implementation checks pass (38 assertions), including independence
of this swap policy to `1e-10` on the four-site untruncated benchmark.
Their output is saved in `results/test-output.txt`.

The results establish a working diagnostic and expose substantial
compression bias. Apparent negative `G` values and small `C` values
must not be promoted to evidence for an all-scale inequality. The
next computational requirement is a contraction method that resolves
the eight-spin reference and retains that accuracy as width increases.

Reproduce the additional comparisons and summary with:

```sh
julia --startup-file=no --project=aklt aklt/scripts/run.jl \
  aklt/configs/thermal-reference-refinement.toml aklt/results/thermal-reference-refinement.json
julia --startup-file=no --project=aklt aklt/scripts/run.jl \
  aklt/configs/thermal-swap-refinement.toml aklt/results/thermal-swap-refinement.json
```

The summary `pilot-summary.json` and the exact projection checks were produced by local Python
helpers (`analyze.py`, `check_projection_*.py`) that are not part of this repository, and the
result files under `aklt/results/` are kept locally too; rerun the configurations above to
regenerate them.
