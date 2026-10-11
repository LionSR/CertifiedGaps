# Staggered-field spin-1/2 Heisenberg certificate

See [REPORT.md](REPORT.md) for results, proof details, conventions and limitations.
See [AUDIT.md](AUDIT.md) for the independent certification audit and saved-artifact checks.
The model is the periodic spin-1/2 Heisenberg chain with staggered Sz field.
An odd-sublattice pi rotation about y gives the exact real RP operators
`{Sx, iSy, Sz}` and the uniform real on-site term `-h Sz`.
For `h=1`, the verified seed `(n0,beta)=(8,21/4)` gives
`gap >= 1.148261490` for **every even L >= 8 sites**. Its final verification run
took about 84 seconds including JSON output. Positive certificates also cover
`h=3/4,1/2,7/20,1/4`. Eigenvalue enclosures, rational seed defects and the recurrence
are rigorous; reference gaps and twenty-site searches are numerical only.

All commands run from the repository root with Julia 1.12.4 and
`OPENBLAS_NUM_THREADS=1`:

```sh
julia --project=julia certificates/heisenberg_staggered_field/scan.jl 16
julia --project=julia certificates/heisenberg_staggered_field/certify.jl '1//1' 16
julia --project=julia certificates/heisenberg_staggered_field/test/runtests.jl
```

`scan.jl 20 '1//4'` extends one field to twenty sites. Blocks larger than 2500
are explored by KrylovKit, explicitly numerical only. `certify.jl` always uses
complete dense eigenvectors in every representative block, with exact symmetry
multiplicities and an asserted dimension sum; it never uses KrylovKit.

The fast checker verifies RP, small-ring spectra against the independent core ED,
sector counts, stored eigenvalue-enclosure propagation and the rational recurrence
inputs. It does not repeat the expensive eigenvector residual calculations; those
are reproduced by `certify.jl`. Stored JSON records include source hashes and versions.

Cloud setup used `/tmp/julia-1.12.4/bin/julia` and
`JULIA_DEPOT_PATH=/tmp/julia-depot`. The official juliaup installer returned HTTP 403;
the official Julia binary archive worked. No existing repository files are changed.
