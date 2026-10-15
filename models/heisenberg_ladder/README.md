# Certified Heisenberg ladder

Spin-1/2 two-leg ladder, leg coupling 1, periodic, one rung per d=4 cell.
The fresh certificate is [results/certified_suite.json](results/certified_suite.json).
[REPORT.md](REPORT.md) gives the proof and checker details; [IDEAS.md](IDEAS.md)
records the routes considered before fresh computation.

**At J_perp=4, Delta >= 1.102403091 for every N=8,12,16,... spin sites.**
Equivalently, every even L>=4 rungs has a unique ground state and this gap.
The seed is (n0,beta0,eps_s,eps_t) =
(4, 2, 3508321//100000000, 123121//5000000).
The survey's rounded 1.102404 is not used as a certified lower bound.

All positive bounds below are freshly certified; all displayed defect decimals
are exact upper bounds, not rounded numerical estimates. Every positive row
uses n0=4 and the same site range N=8,12,16,... .

| J_perp | beta0 | eps_s | eps_t | Certified gap >= |
|---|---:|---:|---:|---:|
| 4 | 2 | 0.03508321 | 0.02462420 | 1.102403091 |
| 7/2 | 9/4 | 0.05806853 | 0.02987064 | 0.778648190 |
| 3 | 5/2 | 0.09896358 | 0.04497730 | 0.454189334 |
| 5/2 | 3 | 0.18558434 | 0.05114736 | 0.042433721 |
| 2 | every beta>0 | No contracting seed at n0=2 or 4 | | |
| 3/2 | every beta>0 | No contracting seed at n0=2 or 4 | | |
| 163069/65536 | 3 | 0.18772796 | 0.05269516 | 0.001531165 |

The negative rows are rigorous **failures of these seeds under the criterion**,
including continuous beta, not claims that the ladder is gapless. Unequal
defects above the diagonal sufficient cutoff are accepted only after the full
core recurrence closes.

## Threshold and ring-size limit

On the grid n0 in {2,4}, beta=k/4 (k=1,...,64), all 128 seeds are rigorously
rejected at J=326133/131072, while a seed is certified at J=163069/65536.
This brackets a local grid transition by the exact endpoints

    2.48819732666015625 < transition < 2.4882354736328125.

This is not a global minimum over continuous beta or J, and monotonicity in J
is not asserted. It does not determine the physical gap-closing coupling.
The lower endpoint is rejected only on the stated beta grid, unlike the
all-beta exclusions at J=2 and 3/2.

The requested historical sector-resolved dense scans through ten rungs
(20 spins, 1,048,576 states) are retained and freshly replayed by
`replay_scans.jl`, including reconstruction of every exact symmetry-sector
multiplicity. Their log Z grids and finite-ring gaps remain floating point;
they were not recomputed in this run. Since n0 must be even, ten-rung traces
do not make n0=5 an eligible seed. Historical heuristic extrapolations suggest
n0=6 at J=2 and n0=8 at J=3/2, requiring 24 and 32 spins respectively; those
estimates are unverified and are not used by any proof.

## Reflection positivity and routes beyond ED

The exact checkerboard rotation is u=exp(i pi Sy)=2iSy on the stated sublattice.
With X=Sx, A=iSy, Z=Sz, all real rational matrices,

    -h_cell = sum_C C tensor C,
    C = {X tensor I, A tensor I, Z tensor I, I tensor X, I tensor A, I tensor Z},
    onsite = -J (X tensor X + A tensor A + Z tensor Z).

`rp.jl` checks both alternating orientations and the original-to-rotated
identities; the suite stores the matrices and also checks them with the core.
This basis change is periodic only for an even number of rungs.

[PERTURBATION.md](PERTURBATION.md) proves that the centered independent-rung
norm sandwich fails for every n and beta when 0<J<=4. This is a documented
negative result for a route reaching arbitrary seed sizes, not an untested
suggestion. Local perturbation is useful: the suite additionally certifies
uniform gap bounds 1.063075969, 0.751905791 and 0.435856503 throughout
|J-J*|<=1/10000 around J*=4,7/2,3 respectively, on the same site range.

## Reproduce

Use Julia 1.12, this model's own project, and one BLAS thread. From the repository root:

```sh
export OPENBLAS_NUM_THREADS=1
julia --project=models/heisenberg_ladder -e 'using Pkg; Pkg.instantiate()'
timeout 2700 julia --project=models/heisenberg_ladder models/heisenberg_ladder/test/runtests.jl
timeout 2700 julia --project=models/heisenberg_ladder models/heisenberg_ladder/certify_suite.jl
timeout 2700 julia --project=models/heisenberg_ladder models/heisenberg_ladder/check.jl
timeout 2700 julia --project=models/heisenberg_ladder models/heisenberg_ladder/replay_scans.jl
timeout 2700 julia --project=models/heisenberg_ladder models/heisenberg_ladder/audit.jl
```

The fresh eight-coupling suite took 366.7 s (6.1 min); its J=4 calculation took
46.1 s. The full fresh checker passed in 308.34 s (5.1 min).
The separate directed-MPFR/trace-moment audit took 37.6 s, and the
historical scan replay took 4.5 s. These recorded times exclude Julia startup
and package installation. All runs have a 2700 s cap. The 57 fast tests pass.
The original ten-rung dense runs took 22.3–33.7 min per coupling; their timings
are historical metadata, not the duration of the current replay.

`check.jl` rejects source/environment hash mismatches, verifies exact RP,
recomputes complete spectra and all seed inequalities, checks the core bound,
replays all 1024 grid decisions, verifies continuous-temperature coverage and
its infinite-beta tail, and checks coupling continuation. No old floating-point
seed or spectrum is trusted. Exact proof numbers use rational strings and
`CertifiedGaps.Results.parse_exact_json`; provenance uses the core result tools.

The proof assumes the repository's abstract criterion and IEEE round-to-nearest
OpenBLAS arithmetic, with explicit outward rounding and error allowances.
Remaining floating point: historical ten-rung scan data, heuristic seed-scale
extrapolations, and timing metadata. Eigensolver output in the fresh certificate
is tentative input to verified enclosures, not an unverified inequality.

`results/certificate.json` and `certify.jl` preserve the older single-J=4 workflow;
`results/certified_suite.json`, `certify_suite.jl`, and `check.jl` are the current
suite. No shared-core or other-model changes are required.
