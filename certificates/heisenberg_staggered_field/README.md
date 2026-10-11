# Staggered-field spin-1/2 Heisenberg certificate

See [REPORT.md](REPORT.md) for results, proof details, conventions and limitations.
See [AUDIT.md](AUDIT.md) for the independent certification audit and saved-artifact checks.
The model is the periodic spin-1/2 Heisenberg chain with staggered Sz field.
An odd-sublattice pi rotation about y gives the exact real RP operators
`{Sx, iSy, Sz}` and the uniform real on-site term `-h Sz`.
For `h=1`, the verified seed `(n0,beta)=(8,21/4)` gives
`gap >= 1.148261` for **every even L >= 8 sites**. Its final verification run
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

For a full spectral replay directly from the saved certificates, without reading
or repeating a parameter search, run:

```sh
julia --project=julia certificates/heisenberg_staggered_field/spectral_replay.jl
```

This reconstructs the exact model blocks, computes fresh complete dense eigenvectors,
and recomputes residuals, orthogonality defects, radii, partition intervals, rational
seeds and bounds for all five certificates. It preserves the originals and writes
`results/spectral_replay.json` with exact comparisons against every stored spectral
field. Binary differences and support of the original nine-decimal claims are
diagnostic and are preserved. Final published bounds are rounded DOWN to six
decimals from the smaller original and fresh outward bound. It is a fresh
execution of the same verifier, not a second independent
implementation of spectral certification.

`spectral_replay_attempt1.jl` and `spectral_replay_attempt2.jl` are archived earlier
scripts with their outputs and failure logs; do not use them as the current replay
command. Attempt 1 imposed binary equality and stopped at a seed mismatch before
saving that field's diagnostics. Attempt 2 preserved fresh data and found the
fresh h=3/4 bound slightly below the old nine-decimal display. The bounded replay
was resumed with `--resume-attempt2`, retaining its two completed fresh spectral
witnesses and computing the remaining three. See AUDIT.md for the exact scope.

Cloud setup used `/tmp/julia-1.12.4/bin/julia` and
`JULIA_DEPOT_PATH=/tmp/julia-depot`. The official juliaup installer returned HTTP 403;
the official Julia binary archive worked. No existing repository files are changed.
