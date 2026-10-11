# Two-leg transverse-field Ising ladder

Certified at **g=3** using **(n0,β0)=(4,3/2)**: **Δ_L ≥ 0.636386863** for every even L≥4 rungs (at least eight spin sites, in multiples of four). This is the smallest successful field in the tested grid. The exact RP basis is σᶻ, with C1=Z⊗I and C2=I⊗Z, each with coefficient one. Verification took about 11 minutes; the 69-test checker passed in 19.93 seconds. The other scans and reference gaps are numerical.


See [REPORT.md](REPORT.md) for the result, assumptions, and scope.

Run from the repository root with Julia 1.12.4 and instantiated `julia/` dependencies:

```sh
export OPENBLAS_NUM_THREADS=1
julia --project=julia certificates/two_leg_tfim/scan.jl 2//1
julia --project=julia certificates/two_leg_tfim/scan.jl 5//2
julia --project=julia certificates/two_leg_tfim/scan.jl 3//1
julia --project=julia certificates/two_leg_tfim/scan.jl 4//1
# certify.jl takes one of the above field values and its stored winning even seed.
julia --project=julia certificates/two_leg_tfim/certify.jl 3//1
julia --project=julia certificates/two_leg_tfim/test/runtests.jl
```

The selected cloud environment has a read-only home directory. Here the official Julia
binary was installed at `/tmp/julia-1.12.4/bin/julia`, and the package depot is
`JULIA_DEPOT_PATH=/tmp/tfim-julia-depot`. Prepend that binary directory to PATH and export
the depot path when reproducing in this environment. The approved juliaup installer
returned HTTP 403; the official versioned Julia binary download succeeded.

No Lanczos or ARPACK is used in the seed partition functions. Scan JSON is numerical;
only `certificate_g*.json` contains verified spectral enclosures. The fast checker
recomputes the exact RP identity and the outward-rounded recurrence from rational seeds;
full residual reproduction requires rerunning `certify.jl`.
