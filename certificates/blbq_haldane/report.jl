using JSON3, Printf
include("result_io.jl")
const ROOT=@__DIR__
d=JSON3.read(read(joinpath(ROOT,"results/scaling.json"),String))
complete14=any(r->r.n==7,d.grid_maxima)
scope=complete14 ? """
Complete dense spectra reach **14 sites** for every requested angle, with numerical
low-energy references through **16 sites**. All 120 (S,k) slots were collected and their
multiplet-weighted dimensions sum to **4,782,969 = 3^14**. This is a numerical negative
grid survey, not a proof that no seed exists at larger eligible scales.
""" : """
Complete dense spectra reach **12 sites**, with numerical references through **16 sites**.
The constructive SU(2) reduction passes small-spectrum tests and 14-site pilots; its
complete 14-site sweep is still being collected. No partial sector set is used as Z14.
"""
f(x)=x===nothing ? "—" : @sprintf("%.6g",x)
open(joinpath(ROOT,"REPORT.md"),"w") do io
    print(io,"""
# Spin-1 bilinear–biquadratic Haldane survey

## Outcome and scope

**No eligible seed was found on the computed grid. No uniform spectral gap is certified.**
$(scope)

One unit cell is one spin-1 site. The tested Hamiltonian is
`H_L = sum_j [cos(theta) X_(j,j+1) + sin(theta) X_(j,j+1)^2]`,
`X = S·S`, with periodic boundaries and angles
`theta = 0, -1/10, -3/10, -1/2, -3/5, -7/10` radians.
For L=2 both directed ring bonds are included, as in the core convention.
There is no on-site term. No energy shift or rescaling is applied.

The scope note's theorem requires **even n0**. Every reported odd n0 is only a
finite-ring diagnostic; odd-n spatial ratios need not be purities of a positive operator.
The hypothetical validity range would be **even L ≥ n0 cells = sites**. Since there is
no certified seed, **there is no asserted validity range or positive bound-to-gap ratio**.

## Exact reflection positivity

Use Cartesian generators `(G_a)_(bc)=epsilon_(abc)`; explicitly

```
G1 = [0 0 0; 0 0 1; 0 -1 0]
G2 = [0 0 -1; 0 0 0; 1 0 0]
G3 = [0 1 0; -1 0 0; 0 0 0]
```

Here `S_a=-i G_a`. Define the **symmetrized**, real symmetric quadrupoles
`Q_ab = -(G_a G_b + G_b G_a)/2` (not the unsymmetrized products).
With `A=sum_a kron(G_a,G_a)` and `B=sum_ab kron(Q_ab,Q_ab)`, the exact identities are

```
X = -A
X^2 = A/2 + B
-h = (cos(theta)-sin(theta)/2) A - sin(theta) B.
```

`model.jl:exact_rp` verifies the two coefficient identities using `Rational{BigInt}`;
this proves the identity for symbolic real c=cos(theta), s=sin(theta), without rounding
transcendental coefficients to rational numbers. For nonzero rational angles these
coefficients are not being claimed to be algebraic. The generator and quadrupole entries
are exact integers/half-integers. The on-site zero matrix is exactly real symmetric.
For -pi/4 < theta ≤ 0, sin(theta) ≤ 0 and cos(theta) > 0, so both weights are nonnegative.

The skew-symmetric and symmetric operator spaces are orthogonal. The three G's span the
former, and six independent Q's span the latter, as checked by a nonzero **exact rational
9-by-9 determinant** in the tests. Consequently the realigned bond is PSD iff the two
weights are nonnegative: a negative weight has a negative direction in its own orthogonal
space. This proves the iff condition for this Cartesian bond, rather than relying on a
floating-point RP test or the numerical converse in the scope note.
An independent exact realignment identity gives eigenvalues `2c-s` (multiplicity 3),
`-s` (multiplicity 5), and `-4s` (multiplicity 1), where c=cos(theta), s=sin(theta).

## Complete spectra and symmetry counting

All partition functions below use every eigenvalue from dense real symmetric LAPACK
solves. Through L12, the reduction uses total Sz=M, cyclic translation k, spin reversal within M=0,
and reflection at k=0 and k=pi. Opposite M and conjugate k share spectra, with explicit
multiplicities 2 where appropriate. At general k, reflection followed by conjugation is
an antiunitary involution; its fixed basis makes the Hamiltonian real. It does not itself
halve the block. Sparse involution projectors supply the additional ± sectors.

The L14 construction instead uses total-spin highest-weight spaces at fixed momentum,
with multiplet weight 2S+1 and conjugate-momentum weight 1 or 2. It does not further split
reflection parity at k=0 or pi. The collector requires all 120 (S,k) entries, including
empty sectors, and checks each dimension against an independent integer character formula
before allowing L14 into a partition function.

Each result records every block dimension and multiplicity, and asserts both
`sum(dim * multiplicity) == 3^L` and that the expanded spectrum has exactly `3^L`
states. SU(2) outputs retain energies with exact integer multiplet/momentum weights.
The tests compare complete reduced and unreduced spectra at L=2,3,4,5 for three
angles (including -0.7). Matrix Hermiticity, real-basis consistency, symmetry invariance
and involution identities are also checked during construction. Complete spectra also
pass analytic first- and second-moment checks at infinite temperature: for L>=3,
`mean(H)=4 L sin(theta)/3` and
`var(H)=L [4 cos(theta)^2/3 - 4 cos(theta)sin(theta)/3 + 8 sin(theta)^2/9]`;
the L=2 variance has an extra factor two because both bonds coincide. These floating checks
validate the numerical implementation; they are not verified eigenvalue enclosures.

| L (sites) | Total dimension | Largest dense block | Stored blocks | Wall seconds, six angles |
|---:|---:|---:|---:|---:|
""")
    for file in sort(filter(f->startswith(f,"spectrum_L"),readdir(joinpath(ROOT,"results")));by=f->parse(Int,match(r"L(\d+)",f)[1]))
        a=read_result_json(joinpath(ROOT,"results",file))
        println(io,"| $(a.L) | $(a.dimension_sum) | $(maximum(s.dimension for s in a.sectors)) | $(length(a.sectors)) | $(f(a.seconds)) |")
    end
    print(io,"""

## Purity scaling table

For each n the entry is the **grid maximum** of min(Ps(n,beta),Pt(2n,beta)), with its
maximizing beta. Grid: beta = 1/8, 2/8, ..., 64 (512 points), using core `scan_seeds`
and `best_seed`. Maxima are numerical, not upper bounds on the continuous-beta optimum.
All even grid points were checked with the actual asymmetric recurrence, not just the
symmetric 0.12731 sufficient threshold. No eligible positive recurrence bound appeared.

| theta | n0 | eligible even seed? | beta at max | max min(Ps,Pt) | Ps | Pt |
|---:|---:|:---:|---:|---:|---:|---:|
""")
    for r in d.grid_maxima
        println(io,"| $(f(r.theta)) | $(r.n) | $(r.eligible ? "yes" : "no") | $(f(r.beta)) | $(f(r.max_min)) | $(f(r.spatial)) | $(f(r.thermal)) |")
    end
    print(io,"""

$(complete14 ? "The n0=7 row uses complete Z14 but is ineligible by parity. No n0=8 grid maximum is reported because full Z16 was not computed." : "No n0=7,8 grid maximum is available until complete Z14/Z16 spectra exist.")
In particular,
the larger-ring Krylov results are **not** silently treated as complete spectra or used
in `logZ_bounds`. Odd-n near-unit diagnostic ratios cannot establish the theorem.

### Best eligible point per parameter (ranked by balanced purity)

These are failed candidate seeds, not certificates. Their stored rational defects are
upward decimal roundings of numerical defects, **not upper enclosures of true defects**.
The checker recomputes `certified_gap_bound` on those exact rational inputs and verifies
that `.bound` is `nothing`. No rational seed is advertised as validated.

| theta | n0 | beta | epsilon_s | epsilon_t | certified bound | bound / ED gap |
|---:|---:|---:|---|---|:---:|:---:|
""")
    for r in d.best_even
        println(io,"| $(f(r.theta)) | $(r.n) | $(f(r.beta)) | $(r.epsilon_s_rational) | $(r.epsilon_t_rational) | none | N/A |")
    end
    evenrows=filter(r->r.eligible,d.grid_maxima)
    closest=argmax(r->r.max_min,evenrows)
    println(io,"\nThe closest computed eligible point is **theta=$(f(closest.theta)), n0=$(closest.n), beta=$(f(closest.beta))**, with balanced purity **$(f(closest.max_min))**. This ranking is numerical and limited to the stated grid.\n")
    print(io,"""
## Ground-energy mismatch and ED reference gaps

`delta(n)=E0(2n)-2 E0(n)` includes the finite-size nonextensive ground-energy contribution.
Positive even-n delta makes Ps fall as exp(-beta delta) in the ground-state-dominated
regime. This competes with the larger beta required for thermal purity. Odd-n deltas are
shown separately through the parity column and must not be fit together with even n.
The gap column uses complete-spectrum E1-E0 where available. Otherwise it is the
numerical M=0 multiplet spacing described below, not a bulk gap assertion.

At L=14 and 16, `reference.jl` searches all momentum and spin-reversal/reflection blocks
in M=0, with six converged low levels per block. Every integer-spin SU(2) multiplet has an
M=0 state, so this covers low multiplet energies numerically. Krylov convergence and
residual norms are recorded, but multiplet ground-state degeneracy is not resolved by
M=0 alone. This spacing equals the full spectral gap if the ground state is a unique
singlet; that condition is not certified by the reference computation. Complete spectra
are preferred when available. These are numerical ED references only, never partition-
function inputs or uniform gap lower bounds.

| theta | n | parity | delta(n) | ED gap at 2n | full spectrum at 2n? |
|---:|---:|:---:|---:|---:|:---:|
""")
    for r in d.energy_scaling
        println(io,"| $(f(r.theta)) | $(r.n) | $(iseven(r.n) ? "even" : "odd") | $(f(r.delta)) | $(f(r.gap_2n)) | $(r.complete_spectrum_2n ? "yes" : "no, reference only") |")
    end
    print(io,"""

## Conditional required-scale forecasts

These are **model-dependent extrapolations, not evidence of an actual seed at those sizes**.
Fit the positive even-n mismatches (n≥4) to both `A exp(-n/xi)` and
`A exp(-n/xi)/sqrt(n)`. Independently extrapolate the last three even finite-size gaps with
`Delta(L)=Delta_inf+a/L` and `Delta(L)=Delta_inf+a/L^2`. The finite-size data are short,
especially near -pi/4, so none of these extrapolations has controlled asymptotic error.
A negative extrapolated gap produces no forecast, rather than being clipped positive.

For an assumed gap g, use the independent-triplet ansatz
`z_L(beta)=(1+3 exp(-beta*g))^L` after removing the ground energy. This gives
`Ps=exp(-beta*delta(n))` and
`Pt=[(1+3 exp(-2 beta*g))/(1+3 exp(-beta*g))^2]^(2n)`.
Scan even n=8:2:2000 and beta*g=2:0.125:30, and report the first n at which the
actual recurrence returns a positive conditional bound. This ansatz ignores dispersion,
interactions, and the detailed low-energy density of states; it is a transparent scale
estimate, not a fit to the full partition function. The variation across four choices is
model sensitivity, **not a confidence interval**.

| theta | gap fit | assumed bulk gap | delta fit | fitted xi | first n0 | beta0 |
|---:|---|---:|---|---:|---:|---:|
""")
    for r in d.extrapolations, e in r.estimates
        for (i,forecast) in enumerate(e.forecasts)
            if forecast===nothing
                println(io,"| $(f(r.theta)) | $(e.method) | $(f(e.gap)) | $(i==1 ? "exponential" : "exponential / sqrt(n)") | — | no forecast | — |")
            else
                println(io,"| $(f(r.theta)) | $(e.method) | $(f(e.gap)) | $(forecast.style) | $(f(forecast.correlation_length)) | $(forecast.n) | $(f(forecast.beta)) |")
            end
        end
    end
    print(io,"""

For theta=-0.5,-0.6,-0.7, the incompatible signs of the two extrapolated gap intercepts
leave the required scale **unresolved**. The finite values shown for those angles are
only outputs of the positive-intercept ansatz. The short series cannot distinguish
the assumed exponential mismatch tail from a long preasymptotic regime.

## Resources and reachable sizes

The environment has Julia **1.12.4**, a **16 GiB** cgroup memory limit and **4 CPU** quota.
The juliaup installer returned HTTP 403. The official Julia 1.12.4 S3 tarball succeeded;
`Pkg.instantiate()` and `Pkg.precompile()` completed using the unmodified core project
and manifest. The runtime/depot are ignored and confined to this certificate directory.
All runs set `BLAS.set_num_threads(1)` and `OPENBLAS_NUM_THREADS=1`; dense parameter solves
use four Julia threads. All changes relative to the base are confined to this certificate directory.

At L=14 the initial M=1,k=1 block has **41,835** states before SU(2) reduction.
One real dense matrix requires **13.04 GiB**; matrix construction plus eigensolver copies
in the implemented path requires more than the available memory (two matrices alone:
**26.08 GiB**). Cubic work is approximately **376 times** the corresponding L=12 block.
At L=16 that block has **310,572** states and needs **718.65 GiB** for one matrix.
`preflight.json` supplies exact integer character counts and the resource projections.

A numerical SU(2) highest-weight construction is now implemented in `su2.jl`.
Its polynomial Casimir projector is checked exactly on small charge spaces, and weighted
complete spectra agree with unreduced dense ED through L=7. A Householder-QR pilot at
L=14, S=2, k=1 produced an 8,106-state dense block in 649.6 seconds including one angle's
eigensolve. The faster CholeskyQR2 version passed the same small tests; its first full
L=14, S=0, k=0 block took 61.8 seconds for all six angles. No incomplete sector collection
is used for Z14. The full sweep is resumable with `run_su2_grid.sh`, with a
45-minute timeout per sector process. See `review_notes.md` and the SU2 pilot JSON.
The exact character counts give largest `(S,k)` highest-weight multiplicities **8,482**
at L=14 (0.536 GiB per dense matrix) and **59,143** at L=16 (26.06 GiB). These counts use
`d(S,k)=d(M=S,k)-d(M=S+1,k)` and assert `sum_(S,k) (2S+1)d(S,k)=3^L`.
Thus this work does **not** establish that full 14-site ED is intrinsically infeasible:
the former basis-construction blocker has been resolved numerically.
$(complete14 ? "The full 14-site sweep is complete. At 16 sites, even the largest SU(2) block's single dense matrix exceeds this environment's memory limit." : "The remaining full sector sweep is the current blocker to the larger-ring purity table.") No
memory-overcommitted 14-site dense solve or incomplete-spectrum purity calculation was run.
With one-site cells, full Z14 adds the odd-n0=7 diagnostic; the next eligible seed
beyond n0=6 is n0=8, which needs complete Z16 under the stated two-purity workflow.

## Rigor boundary and eigenvalue-enclosure theorem

**Proved exactly:** the RP coefficient identities, real on-site term, and coefficient
signs in the stated parameter window. Integer Hilbert-space counts are asserted.
**Numerical only:** matrix symmetry checks, every computed spectrum, all purities, failed
seed defects, finite-size gaps, delta fits, rankings, and extrapolated required scales.
**Not established:** a true-spectrum eigenvalue enclosure, an admissible rational seed,
a uniform gap, uniqueness from the criterion, or any thermodynamic extrapolation error bar.

The requested best-seed rigorous rerun is not triggered because no eligible candidate
passes. Consequently no residual/orthogonality eigenvalue enclosure is claimed or
implemented. For clarity, the theorem such a future path must apply is the polar-factor
Weyl bound: for square approximate eigenvectors Q and real diagonal Lambda, if verified
`||Q'*Q-I||_2 <= eta < 1`, `||H*Q-Q*Lambda||_2 <= rho`, then there is a unitary polar factor U
with `||U-Q||_2 <= 1-sqrt(1-eta)`, and sorted eigenvalues of H lie within
`rho+(||H||_2+||Lambda||_2)*(1-sqrt(1-eta))` of sorted Lambda. This follows by bounding
`U'*H*U-Lambda` and applying Weyl. All norms, matrix assembly errors, trigonometric
coefficients and arithmetic errors must be enclosed outward. Merely accumulating in
BigFloat without a rounding bound is not sufficient.

Given eigenvalue intervals, logZ lower and upper bounds follow by summing exp(-beta E)
with directed rounding, and purity lower bounds use
`exp(logZ2n_lo-2 logZn_hi)` and `exp(logZ2n_2beta_lo-2 logZ2n_beta_hi)`.
Only then may rational upper defect bounds be passed to the certified recurrence.
The numerical `PartialSpectrum.exact=true` flag here means *complete*, not *certified*.

## Reproduction and validation

From the repository root (the wrapper selects the downloaded Julia and private depot):

```sh
certificates/blbq_haldane/run.sh -t 4 certificates/blbq_haldane/scan.jl 2 3 4 5 6 7 8 10 12
certificates/blbq_haldane/run.sh certificates/blbq_haldane/reference.jl 14
certificates/blbq_haldane/run.sh certificates/blbq_haldane/reference.jl 16
certificates/blbq_haldane/run.sh certificates/blbq_haldane/preflight.jl
certificates/blbq_haldane/run.sh certificates/blbq_haldane/analyze.jl
certificates/blbq_haldane/run.sh certificates/blbq_haldane/report.jl
certificates/blbq_haldane/run.sh certificates/blbq_haldane/test/runtests.jl
```

Each invocation runs `julia --project=julia`; each completed numerical run is below the
45-minute cap. Stored timers cover the calculation functions, excluding environment setup
and initial package imports. Result JSON uses the shared core's `Results.write_result`
(ThermalMPOGap historically, CertifiedGaps after the rename). Large complete spectra are
losslessly gzip-compressed; `archive.json` records original and compressed hashes.
Result metadata contains source hashes and run metadata. Exact rational tokens are strings. Spectrum
files retain every eigenvalue and every sector dimension, enabling independent numerical
reanalysis without repeating ED. `scaling.json` hashes its spectrum/reference inputs.
The fast test recomputes exact RP, small-system spectrum crosschecks, complete dimension
counts, and the failed recurrence outputs from stored exact rational tokens. It does not
pretend to certify the rejected numerical input defects. See `validation_current_main.json`
for current-main test timing and `historical/thermal_validation.json` for the original
snapshot. `MIGRATION.md` describes the independent source/dependency comparison. No merge
has been performed; draft publication is tracked in `PR.md` when available.
""")
    print(io,"\n### Larger-ring reference runtimes\n\n| L | Wall seconds |\n|---:|---:|\n")
    for L in (14,16)
        a=JSON3.read(read(joinpath(ROOT,"results/reference_L$(L).json"),String))
        println(io,"| $L | $(f(a.seconds)) |")
    end
end
