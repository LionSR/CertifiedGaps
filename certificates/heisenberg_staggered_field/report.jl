include("common.jl")
using Printf
fmt(x)=@sprintf("%.9g",x)
function boundfmt(x)
    a=floor(BigInt,Rational{BigInt}(x)*big(10)^6)
    string(a÷10^6,".",lpad(string(a%10^6),6,'0'))
end
readjson(path)=JSON3.read(read(path,String),Dict{String,Any})
certs=[readjson(joinpath(ROOT,"results","certificate_$(tag(h)).json")) for h in FIELDS]
replay=readjson(joinpath(ROOT,"results","spectral_replay.json"))
@assert replay["status"]=="complete"
for c in certs
    fresh=only(filter(r->r["h"]==c["h"],replay["records"]))
    c["bound"]=min(c["bound"],fresh["bound"])
    c["ratio"]=c["bound"]/c["reference_gap"]
end
scans=[readjson(joinpath(ROOT,"results","scan_$(tag(h))_20.json")) for h in FIELDS]
threshold=[readjson(joinpath(ROOT,"results","scan_$(tag(h))_20.json")) for h in [1//5,3//20,1//10]]
open(joinpath(ROOT,"REPORT.md"),"w") do io
    print(io,"""
# Uniform gap certificates: Heisenberg chain in a staggered field

## Outcome and scope

All five requested fields have rigorous certificates from 8- and 16-site rings.
The final bounds apply to **every even L >= 8 one-spin cells = sites**. The physical
Hamiltonian has a two-site field pattern, but its rotated representative has a
one-spin unit cell. The theorem only uses even n0; odd n0 are not valid seeds here.
No threshold occurs within the requested field set: the smallest certified field
is 1/4. This is a certificate/search threshold statement, not a physical critical field.

A separate numerical search reaches 20 sites using Sz and translation-by-two
momentum sectors. It finds stronger candidates but uses partial Krylov spectra
for blocks above dimension 2500; these candidates are **numerical, not certified**.
We certify the best seed on the complete dense grid through 16 sites, rather than
claim certification of the best candidate on the larger, partial-spectrum grid.
The requested fully dense certification of the best 20-site candidates remains
unperformed. No 24-site computation is claimed.
The largest block grows from 1620 to 18504: a single complex dense matrix then
occupies about 5.48 GB, and the current verifier retains several such arrays.
The naive cubic eigensolve cost grows by about 1490 times per largest block.
These are resource estimates, not measured failed dense twenty-site runs.

## Model, exact RP identity, and length convention

H = sum_j S_j dot S_(j+1) - h sum_j (-1)^j Sz_j, with periodic bonds,
spin 1/2, exchange strength 1, and h in {1,3/4,1/2,7/20,1/4}.
The two-site ring includes both periodic bonds, as does the core implementation.
Define exact rational matrices

```
X = [0  1/2;  1/2  0],
J = [0  1/2; -1/2  0] = i Sy,
Z = [1/2 0; 0 -1/2],   U = 2J.
```

Rotate odd sites by the real orthogonal U (a pi rotation about y). Both orientations
of an adjacent rotated/unrotated pair give

```
h_rot = -X tensor X - J tensor J - Z tensor Z,
-h_rot = X tensor X + J tensor J + Z tensor Z.
```

Every coefficient c is exactly 1 and every C is real. J is antisymmetric; the
hypothesis requires real C, not symmetric C. U Z U^T = -Z, so the staggered field
becomes the uniform real symmetric on-site term -h Z. `exact_rp` checks the bond
identities in Rational arithmetic in both orientations, and checks the on-site
term. Thus the stated spatial transfer self-adjointness hypothesis holds. The
positive two-purity bound then gives uniqueness as well as the gap, using the
repository criterion; no numerical assumption of uniqueness is substituted.

## Dense-grid certificates

Grid: n0 in {2,4,6,8}, beta in {0.25,0.50,...,32}. `scan_seeds` and `best_seed`
select the largest bound. Only that seed is rerun with dense eigenvectors and
outward error bounds for each requested field. All exact rational parameters
and seed defects are strings in the certificate JSON.

| h | n0 | beta | rigorous bound | numerical 16-site gap | bound/gap | certification seconds |
|---|---:|---:|---:|---:|---:|---:|
""")
    for c in certs
        println(io,"| $(c["h"]) | $(c["n0"]) | $(c["beta"]) | $(boundfmt(c["bound"])) | $(fmt(c["reference_gap"])) | $(fmt(c["ratio"])) | $(fmt(c["seconds"])) |")
    end
    println(io,"\nDisplayed rigorous bounds are rounded DOWN to six decimal places from the smaller original and freshly replayed bound. Original raw outputs remain in certificate JSONs; fresh raw outputs and exact seed differences are in `results/spectral_replay.json`. The full spectral replay found small numerical differences, so the earlier nine-decimal display is superseded; see AUDIT.md.")
    c=certs[1]
    println(io,"\nFor h=1, the exact stored upper seed defects are epsilon_s = `$(c["epsilon_s"])` and epsilon_t = `$(c["epsilon_t"])`, with beta = `$(c["beta"])`.")
    println(io,"\n| h=1 ring | max entry error | max residual upper | max orthogonality upper | max eigenvalue radius |")
    println(io,"|---:|---:|---:|---:|---:|")
    for (key,L) in (("small_sectors",8),("large_sectors",16))
        rs=c[key]
        println(io,"| $L | $(join([fmt(maximum(r[k] for r in rs)) for k in ("entry_error","residual_upper","orthogonality_upper","radius")]," | ")) |")
    end
    print(io,"""

The h=1 survey value supplied in the task is 1.147821, with reference gap 1.87351.
The current repository recurrence gives 1.1490006457432576 from unverified dense
seed values at (8,5.25). Our verified bound includes enclosure loss and is the
value above. We do not silently substitute the survey number. Reference gaps and
ratios in this report are numerical ED values, not independent uniform gap proofs.

## Complete sectors and multiplicity

In the physical frame, p counts down spins, m=L/2, and k indexes momentum under
translation by two sites. An orbit of period ell contributes iff k*ell is divisible
by m. Its normalized state is sum_r exp(-2 pi i k r/m) |T^(2r) a>/sqrt(ell).
For a transition from representative a into T^(2t)b, the matrix element is the
local amplitude times sqrt(ell_a/ell_b) exp(2 pi i k t/m). This formula is used
directly; a small-ring test compares the full merged spectrum with the independent
core Sz implementation.

Global spin flip composed with one-site translation maps p to L-p and preserves
H and momentum. Complex conjugation maps k to m-k. We diagonalize one dense
representative per such exact equivalence, with multiplicities 1,2, or 4; no
uncomputed inequivalent sector is omitted. Exact integer orbit counts and the
assertion sum(dimension*multiplicity)==2^L run for every ring. No Krylov solve is
called on the certificate path.

| ring | representative blocks | largest block | weighted dimension | required dimension |
|---:|---:|---:|---:|---:|
""")
    for (key,L) in (("small_sectors",8),("large_sectors",16))
        rs=certs[1][key]
        println(io,"| $L | $(length(rs)) | $(maximum(r["dimension"] for r in rs)) | $(sum(r["dimension"]*r["multiplicity"] for r in rs)) | $(2^L) |")
    end
    println(io,"\nFull sector dimensions, multiplicities, eigenvalues, radii and timings are stored per block in each certificate JSON. For L=16:")
    println(io,"\n| p | dimensions for k=0,1,2,3,4 | charge multiplicity |")
    println(io,"|---:|---|---:|")
    for p in 0:8
        rs=filter(r->r["p"]==p,certs[1]["large_sectors"])
        println(io,"| $p | $(join([string(r["k"])*":"*string(r["dimension"]) for r in rs],", ")) | $(p==8 ? 1 : 2) |")
    end
    print(io,"""

Momentum multiplicities for k=0,1,2,3,4 are 1,2,2,2,1. Absent entries have dimension zero.

## Eigenvalue-enclosure theorem and implementation

All stored Float64 entries are interpreted as exact binary numbers. Let H be the
exact Hermitian momentum block and A the Hermitian matrix submitted to LAPACK.
`block(...; verified=true)` evaluates rational coefficients, orbit normalization
square roots, and trigonometric phases in outward-rounded IntervalArithmetic.
Pi itself is enclosed, not replaced with a single floating-point approximation.
Each exact upper-triangular entry is enclosed; the lower triangle is its conjugate.
This gives a Frobenius upper bound a >= ||H-A||_2, stored as `entry_error`.

Let Q and Lambda be the full computed eigenvector/eigenvalue arrays. The program
explicitly calculates Rhat=fl(AQ-Q Lambda) and Ehat=fl(Q* Q-I). It never assumes
that LAPACK's answer is exact. Each Frobenius norm is bounded using upward scalar
multiplication, addition, and square root (`nextfloat` after each operation,
including a subnormal allowance). Complex real and imaginary components are
accounted for separately. Write q,h,l for these upper norms of Q,A,Lambda.
With u=2^-53 and g=gamma_(32d+32)=(32d+32)u/(1-(32d+32)u), the implementation uses

```
r = ||Rhat||_F,up + g (h q + q l) + 10^-280,
e = ||Ehat||_F,up + g (q^2 + sqrt(d)) + 10^-280.
```

These follow from the standard dot-product rounding bound
|fl(sum products)-sum products| <= gamma_s sum |products|, applied to real and
imaginary components and their subtractions, and || |A| |B| ||_F <= ||A||_F||B||_F.
The factor 32d+32 conservatively exceeds the number of elementary operations per
complex dot product and subtraction. Sparse AQ has no longer dot products than
the dense dimension. The sqrt(d) term explicitly budgets subtraction of I.
The added 10^-280 exceeds the aggregate possible subnormal
absolute errors at all reached dimensions. Inputs/results are finite and the
chosen gamma denominator is positive. No fast-math transformations are requested.
This is a rigorous norm error bound on ordinary matrix multiplication, rather
than merely trusting a residual accumulated in the same precision.

**Polar/Weyl argument (with multiplicity):** assert e<1. The polar unitary V of Q
exists and ||Q-V||_2 <= dQ=e/(1+sqrt(1-e)). Therefore

```
||A - V Lambda V*||_2
 = ||AV - V Lambda||_2
 <= r + (||A||_2 + ||Lambda||_2) dQ.
```

Use ||A||_2 <= L(3/4+h/2)+a and ||Lambda||_2=max|lambda_i|. Weyl's inequality
for the two Hermitian matrices H and V Lambda V* encloses each ordered eigenvalue,
with multiplicity, in lambda_i +/- eta, where

```
eta = a + r + (L(3/4+h/2)+a+max|lambda_i|) dQ.
```

All scalar operations defining r,e,dQ,eta are outward interval operations. This
argument does not rely on simple eigenvalues or on matching individual residuals
to distinct eigenvalues. Every eigenvalue of each block is present.

For beta>0, `enclosed_logz` sums exp(-beta*[lambda_i-eta,lambda_i+eta]) with exact
symmetry multiplicities and a common exact binary shift, then takes an interval
log. Both lower and upper logZ endpoints are retained. The lower purity bounds
come from logZ(2n,beta)-2logZ(n,beta) and logZ(2n,2beta)-2logZ(2n,beta).
Upper endpoints of 1-exp(these intervals) are converted to exact Rational{BigInt}
values, with no intermediate decimal rounding. Those rational seeds and rational
beta are passed to the repository's `certified_gap_bound`, including its rigorous
tail correction. The fast test recomputes this propagation from stored enclosures.

## Twenty-site exploration and threshold

The following results are **numerical, not certified**. Full dense spectra are used
through 16 sites. At 20 sites only blocks of dimension <=2500 are dense; larger
blocks retain 80 low levels from KrylovKit. All sectors are visited, with the same
exact dimension count. A common truncation cutoff and a dimension-weighted tail
are included via core `PartialSpectrum`/`logZ_bounds`, but Krylov convergence and
multiplicities are not verified, so these are not rigorous partition bounds.
Exact twenty-site dimensions are recorded below.

| h | best n0 | beta | numerical candidate bound | numerical 20-site reference gap |
|---|---:|---:|---:|---:|
""")
    for s in vcat(scans,threshold)
        b=s["best"]
        println(io,b===nothing ? "| $(s["h"]) | — | — | no grid seed | $(fmt(s["reference_gap"])) |" :
            "| $(s["h"]) | $(b["n0"]) | $(fmt(b["beta"])) | $(fmt(b["bound"])) | $(fmt(s["reference_gap"])) |")
    end
    r20=readjson(joinpath(ROOT,"results","spectrum_1_4_20.json"))["sectors"]
    println(io,"\nActual L=20 largest block: $(maximum(r["dimension"] for r in r20)); weighted total: $(sum(r["dimension"]*r["multiplicity"] for r in r20)) = 2^20.")
    println(io,"\nAdditional fields 1/5,3/20,1/10 probe the numerical search threshold beyond the requested set. A failed grid point is not evidence of a gapless phase.")
    successes=filter(s->s["best"]!==nothing,vcat(scans,threshold))
    failures=filter(s->s["best"]===nothing,vcat(scans,threshold))
    if !isempty(failures)
        lo=maximum(Float64(rat(s["h"])) for s in failures)
        hi=minimum(Float64(rat(s["h"])) for s in successes)
        println(io,"\nThe sampled twenty-site search fails at h=$lo and succeeds at h=$hi. This brackets an apparent finite-grid threshold between those sampled fields; no monotonicity between them or physical critical field is proved.")
    end
    println(io,"\nFor every explored field, the following gives max_beta min(Ps,Pt) and delta(n)=E0(2n)-2E0(n). Partial-spectrum values remain numerical.")
    println(io,"\n| h | n0 | max min purity | delta(n) |")
    println(io,"|---|---:|---:|---:|")
    for s in vcat(scans,threshold),r in s["per_n"]
        println(io,"| $(s["h"]) | $(r["n0"]) | $(fmt(r["max_min_purity"])) | $(fmt(r["delta"])) |")
    end
    print(io,"""

## Predicted scale, extrapolation, and limits

Using the supplied Delta proportional to h^(2/3) scaling and calibrating the
**quality of the h=1 survey seed**, the naive comparable-quality scales are
n0 approximately 8 h^(-2/3), beta approximately 5.25 h^(-2/3). Logarithmic corrections
and changes in amplitudes are omitted. This predicts a scale, not a necessary
condition for any positive bound. Weaker positive certificates can occur earlier.

| h | predicted n0 (continuous) | next even n0 | predicted beta | largest searched n0 |
|---|---:|---:|---:|---:|
""")
    for h in FIELDS
        f=Float64(h)^(-2/3)
        println(io,"| $h | $(fmt(8f)) | $(2ceil(Int,4f)) | $(fmt(5.25f)) | 10 |")
    end
    for s in threshold
        s["best"]===nothing || continue
        rows=s["per_n"];xs=Float64[r["n0"] for r in rows[end-2:end]]
        ys=log.(Float64[r["delta"] for r in rows[end-2:end]])
        intercept,slope=hcat(ones(3),xs)\ys
        gap=s["reference_gap"];found=nothing
        # Explicit heuristic only: neglect thermal spatial corrections, use 4n independent
        # low-lying excitations at the observed gap. Extrapolation is not a bound.
        for n in 12:2:128
            candidates=scan_seeds((n,b)->(1-exp(-b*exp(intercept+slope*n)),
                1-(1+4n*exp(-2b*gap))/(1+4n*exp(-b*gap))^2),[n],0.25:0.25:128)
            bs=best_seed(candidates)
            if bs!==nothing;found=bs;break;end
        end
        println(io,"\nFor h=$(s["h"]), fitting the last three delta(n) values gives log delta approximately $(fmt(intercept)) + ($(fmt(slope))) n, or correlation scale $(fmt(-1/slope)).")
        if found!==nothing
            println(io,"A heuristic model Ps=exp[-beta delta(n)] and Pt=(1+4n exp[-2beta gap])/(1+4n exp[-beta gap])^2 first yields a predicted seed at (n0,beta)=($(found.n₀),$(found.β₀)), requiring $(2found.n₀) sites. This is an extrapolation, not an existence proof or certificate.")
        else
            println(io,"That heuristic does not yield a seed through n0=128, beta=128.")
        end
    end
    print(io,"""

## Reproduction, environment, and evidence status

Julia 1.12.4 was verified. The approved juliaup installer returned HTTP 403; the
official julialang-s3 binary archive succeeded. Dependencies were instantiated and
precompiled from the unchanged julia project. OPENBLAS_NUM_THREADS=1 and
BLAS.set_num_threads(1) were used. The environment has approximately 17 GiB RAM
and five CPUs. Independent parameter runs used separate processes. No individual
run approached the 45-minute limit. Source SHA-256 hashes, package versions,
per-block times, and Julia/BLAS thread counts are recorded in JSON.

Run commands are in README.md. `test/runtests.jl` reruns the exact RP identities,
compares small complete momentum spectra with the independent core ED, checks
sector totals, recomputes enclosed partition functions, reconstructs rational
seeds, and reruns the certified recurrence. It does not reconstruct the large
Q matrices; rerun `certify.jl` to audit residual calculations from scratch.

**Certified:** exact RP at the five requested rational fields; complete weighted
sector coverage on certificate rings; eigenvalue, partition-function and purity
enclosures under IEEE arithmetic and the stated error analysis; the reported
uniform lower bounds for every even L>=8 sites, with uniqueness from the criterion.

**Numerical only:** grid optimization, all reference gaps and ratios, twenty-site
partial spectra and seed candidates, finite-grid threshold observations,
finite-size delta fits, and h^(2/3) scale predictions.

**Remaining limit:** the stronger best candidates on the twenty-site search grid
have not been rerun through complete dense verification. No physical threshold
or impossibility of certificates at weaker fields is proved. No existing core,
other certificate directory, or notes were edited; nothing was merged or published.
""")
    println(io,"\nActual checker result: **71 tests passed**; $(strip(read(joinpath(ROOT,"results","test_timing.txt"),String))). See `results/tests.log`. Source-hash audit and full dependency versions are in `results/environment.json`.")
    println(io,"\nThe separate independent certification review and its acceptance status are recorded in [AUDIT.md](AUDIT.md).")
end
