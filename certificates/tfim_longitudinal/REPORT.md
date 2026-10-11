# Longitudinal-field TFIM: uniform spectral-gap certificates

**30 of 30 parameter points certified.** One cell is one spin-1/2 site.
Every bound below holds for **every even L ≥ n₀ cells = sites**, not merely the ED sizes.
The grid is g ∈ {0.25,0.5,0.75,1,1.5,2}, h ∈ {0.05,0.1,0.25,0.5,1}.
No result is asserted for an unsampled parameter or for odd L.

## Model, RP and uniqueness

H_L = −Σ Z_j Z_(j+1) − gΣ X_j − hΣ Z_j, with periodic bonds, including both bonds at L=2.
In the σᶻ basis, X=[[0,1],[1,0]], Z=diag(1,−1), c=1 and C=Z.
Thus **−h_bond = 1·Z⊗Z**, with exact rational entries and positive coefficient.
The on-site matrix [[−h,−g],[−g,h]] is real symmetric. `rp_check` verifies these identities
in Rational{BigInt} for every stated rational parameter. The spatial transfer operator is therefore
self-adjoint under the repository's Appendix E theorem. The ground state at every finite L is
unique independently: for g>0, a sufficiently large scalar shift of −H is an irreducible
nonnegative matrix on the connected spin-flip hypercube. Perron–Frobenius gives a simple largest
eigenvalue. That finite-size argument alone does not prove a uniform gap.

## Search and actual sizes

The initial scan uses n₀ ∈ {2,4,6}, β₀ ∈ {0.25,0.50,…,64}, via `scan_seeds` and `best_seed`.
Only even n₀ are eligible: the stated criterion requires an even spatial power.
All lengths in each JSON's `spectra` were densely diagonalized in **every** momentum sector.
No Lanczos, ARPACK, DMRG, or partial-spectrum tail is used. The best numerical seed per point
was rerun with rigorous spectral and partition-function enclosures. A numerical optimum need not
maximize the final conservative certified bound. All 30 points succeeded at this initial size cap; no larger rings were needed for the uniform claims.

| L (sites = cells) | Momentum sector dimensions k=0,…,L−1 | Sum |
|---:|---|---:|
| 2 | 3, 1 | 4 = 2^2 |
| 4 | 6, 3, 4, 3 | 16 = 2^4 |
| 6 | 14, 9, 11, 10, 11, 9 | 64 = 2^6 |
| 8 | 36, 30, 33, 30, 34, 30, 33, 30 | 256 = 2^8 |
| 12 | 352, 335, 344, 338, 346, 335, 348, 335, 346, 338, 344, 335 | 4096 = 2^12 |

The requested 20–24-site full-spectrum survey was **not reached**. At L=20 even generic k=1
has 52,377 states: a single complex dense matrix takes about 40.9 GiB, before eigenvectors and
workspace, on a machine with a 16 GiB container memory limit. At L=24 this rises to thousands of GiB.
See `results/audit.json` for exact integer resource counts. Even a real representation does not
make the 20-site generic block fit with eigensolver workspace. Spin-flip parity is broken by h>0.
Reflection exchanges k and −k, rather than independently halving every block; the present code
does not further split k=0,π by reflection parity. Smaller seeds already establish uniform results.

## Best seed and bound per parameter

Decimals are approximate displays; exact input rationals and exact binary lower bounds are in JSON.
Δ_ED is the numerical first excitation gap at the largest scanned ring, **not** an extrapolated
thermodynamic gap and not a rigorous reference. Ratios compare the certified bound to this ED value.

| g | h | n₀ | β₀ | certified bound ≈ | ED gap ≈ | bound/ED | total seconds |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.25 | 0.05 | 6 | 24.25 | 0.392289 | 1.190350 | 0.329558 | 19.138308 |
| 0.25 | 0.1 | 6 | 13.75 | 0.733243 | 2.380653 | 0.308001 | 8.659905 |
| 0.25 | 0.25 | 6 | 6.5 | 1.660603 | 4.190218 | 0.396305 | 8.698939 |
| 0.25 | 0.5 | 6 | 5.5 | 1.975681 | 4.842757 | 0.407966 | 8.662684 |
| 0.25 | 1.0 | 6 | 5.25 | 2.033620 | 5.939665 | 0.342379 | 8.749649 |
| 0.5 | 0.05 | 6 | 13.0 | 0.456363 | 1.157530 | 0.394256 | 8.741896 |
| 0.5 | 0.1 | 6 | 8.0 | 0.959161 | 2.314523 | 0.414410 | 8.657309 |
| 0.5 | 0.25 | 6 | 4.5 | 2.361602 | 3.641722 | 0.648485 | 8.682797 |
| 0.5 | 0.5 | 6 | 4.0 | 2.775375 | 4.491823 | 0.617873 | 8.687641 |
| 0.5 | 1.0 | 6 | 4.0 | 2.727630 | 5.780397 | 0.471876 | 8.673671 |
| 0.75 | 0.05 | 6 | 8.0 | 0.317029 | 1.080717 | 0.293350 | 8.598307 |
| 0.75 | 0.1 | 6 | 5.25 | 0.771217 | 2.120935 | 0.363621 | 8.474545 |
| 0.75 | 0.25 | 6 | 3.5 | 2.101814 | 3.075642 | 0.683374 | 8.533159 |
| 0.75 | 0.5 | 6 | 3.5 | 3.011730 | 4.089102 | 0.736526 | 8.610338 |
| 0.75 | 1.0 | 6 | 3.5 | 3.137425 | 5.567299 | 0.563545 | 8.697487 |
| 1.0 | 0.05 | 6 | 6.5 | 0.076325 | 0.867634 | 0.087969 | 8.679040 |
| 1.0 | 0.1 | 6 | 4.25 | 0.504270 | 1.532388 | 0.329075 | 8.681801 |
| 1.0 | 0.25 | 6 | 3.0 | 1.609224 | 2.572964 | 0.625436 | 8.659700 |
| 1.0 | 0.5 | 6 | 3.0 | 2.752299 | 3.708922 | 0.742075 | 8.460910 |
| 1.0 | 1.0 | 6 | 3.0 | 3.640306 | 5.342466 | 0.681391 | 8.836568 |
| 1.5 | 0.05 | 6 | 3.5 | 0.169162 | 1.085259 | 0.155873 | 8.360470 |
| 1.5 | 0.1 | 6 | 3.25 | 0.360094 | 1.288131 | 0.279548 | 8.415947 |
| 1.5 | 0.25 | 6 | 2.75 | 1.055623 | 2.062622 | 0.511787 | 9.017051 |
| 1.5 | 0.5 | 6 | 2.5 | 2.129274 | 3.204046 | 0.664558 | 8.348947 |
| 1.5 | 1.0 | 6 | 2.5 | 3.704972 | 4.978175 | 0.744243 | 8.371029 |
| 2.0 | 0.05 | 6 | 2.25 | 0.912083 | 2.019622 | 0.451611 | 8.376775 |
| 2.0 | 0.1 | 6 | 2.25 | 0.964253 | 2.076256 | 0.464419 | 8.389328 |
| 2.0 | 0.25 | 6 | 2.25 | 1.272783 | 2.414984 | 0.527036 | 8.447028 |
| 2.0 | 0.5 | 6 | 2.25 | 2.000505 | 3.222034 | 0.620883 | 8.400095 |
| 2.0 | 1.0 | 6 | 2.0 | 3.491170 | 4.853869 | 0.719255 | 8.609607 |

The h>0 points at **g=1** and in the **ordered zero-field regime g<1** are included explicitly.
The phrase 'critical point' here means the g=1 location of the h=0 model; h>0 is a different,
non-integrable Hamiltonian. Nothing here certifies a uniform gap as h tends to zero.

The unequal-defect recurrence can succeed even when max min(Pₛ,Pτ) is below 1−0.12731; the equal-defect square is sufficient, not necessary.

## Region map

Entries are certified lower bounds, approximately displayed; row g, column h.

| g \ h | 0.05 | 0.1 | 0.25 | 0.5 | 1.0 |
|---:|---:|---:|---:|---:|---:|
| 0.25 | 0.392289 | 0.733243 | 1.660603 | 1.975681 | 2.033620 |
| 0.5 | 0.456363 | 0.959161 | 2.361602 | 2.775375 | 2.727630 |
| 0.75 | 0.317029 | 0.771217 | 2.101814 | 3.011730 | 3.137425 |
| 1.0 | 0.076325 | 0.504270 | 1.609224 | 2.752299 | 3.640306 |
| 1.5 | 0.169162 | 0.360094 | 1.055623 | 2.129274 | 3.704972 |
| 2.0 | 0.912083 | 0.964253 | 1.272783 | 2.000505 | 3.491170 |

## Numerical finite-size diagnostics

Each pair is max over the β grid of min(Pₛ,Pτ), followed by δ(n)=E₀(2n)−2E₀(n).
These diagnostics use unverified numerical eigenvalues; they are not certificate inputs.
| g | h | n : (max min purity; δ) |
|---:|---:|---|
| 0.25 | 0.05 | 2: (0.823887; 8.21818e-03) · 4: (0.999845; 4.71886e-06) · 6: (1.000000; 5.14624e-09) |
| 0.25 | 0.1 | 2: (0.942730; 3.71832e-03) · 4: (0.999970; 1.58777e-06) · 6: (1.000000; 1.22463e-09) |
| 0.25 | 0.25 | 2: (0.990851; 1.07164e-03) · 4: (0.999998; 2.16942e-07) · 6: (1.000000; 6.99085e-11) |
| 0.25 | 0.5 | 2: (0.998289; 3.26691e-04) · 4: (1.000000; 2.38111e-08) · 6: (1.000000; 2.47269e-12) |
| 0.25 | 1.0 | 2: (0.999772; 7.12615e-05) · 4: (1.000000; 1.07761e-09) · 6: (1.000000; 3.19744e-14) |
| 0.5 | 0.05 | 2: (0.515966; 9.41834e-02) · 4: (0.979492; 1.05754e-03) · 6: (0.999629; 1.77961e-05) |
| 0.5 | 0.1 | 2: (0.657471; 5.04400e-02) · 4: (0.995702; 3.54781e-04) · 6: (0.999949; 4.17069e-06) |
| 0.5 | 0.25 | 2: (0.912333; 1.57176e-02) · 4: (0.999691; 4.84429e-05) · 6: (0.999999; 2.35002e-07) |
| 0.5 | 0.5 | 2: (0.980912; 4.90042e-03) · 4: (0.999979; 5.38490e-06) · 6: (1.000000; 8.47346e-09) |
| 0.5 | 1.0 | 2: (0.997226; 1.08926e-03) · 4: (0.999999; 2.50896e-07) · 6: (1.000000; 7.67386e-11) |
| 0.75 | 0.05 | 2: (0.511618; 2.57384e-01) · 4: (0.792806; 1.98354e-02) · 6: (0.977961; 1.62239e-03) |
| 0.75 | 0.1 | 2: (0.546211; 1.73549e-01) · 4: (0.942818; 6.90249e-03) · 6: (0.996667; 3.75230e-04) |
| 0.75 | 0.25 | 2: (0.761125; 6.60480e-02) · 4: (0.995122; 9.66219e-04) · 6: (0.999893; 2.11994e-05) |
| 0.75 | 0.5 | 2: (0.932607; 2.20140e-02) · 4: (0.999625; 1.11181e-04) · 6: (0.999997; 7.96749e-07) |
| 0.75 | 1.0 | 2: (0.988842; 5.09927e-03) · 4: (0.999986; 5.47924e-06) · 6: (1.000000; 7.81831e-09) |
| 1.0 | 0.05 | 2: (0.507059; 3.88036e-01) · 4: (0.604279; 9.30832e-02) · 6: (0.816466; 2.31877e-02) |
| 1.0 | 0.1 | 2: (0.550969; 3.12475e-01) · 4: (0.791579; 3.98586e-02) · 6: (0.960197; 5.92191e-03) |
| 1.0 | 0.25 | 2: (0.648078; 1.53158e-01) · 4: (0.973151; 6.42237e-03) · 6: (0.998359; 3.71077e-04) |
| 1.0 | 0.5 | 2: (0.861184; 5.77726e-02) · 4: (0.997576; 8.04603e-04) · 6: (0.999948; 1.55371e-05) |
| 1.0 | 1.0 | 2: (0.972366; 1.43935e-02) · 4: (0.999900; 4.35900e-05) · 6: (0.999999; 1.74852e-07) |
| 1.5 | 0.05 | 2: (0.578412; 4.42362e-01) · 4: (0.746042; 1.23891e-01) · 6: (0.869695; 4.18772e-02) |
| 1.5 | 0.1 | 2: (0.596205; 4.18526e-01) · 4: (0.779107; 9.89615e-02) · 6: (0.907139; 2.78620e-02) |
| 1.5 | 0.25 | 2: (0.657538; 3.09299e-01) · 4: (0.901249; 3.70886e-02) · 6: (0.981824; 5.43265e-03) |
| 1.5 | 0.5 | 2: (0.762014; 1.63980e-01) · 4: (0.982105; 7.41331e-03) · 6: (0.998752; 4.33227e-04) |
| 1.5 | 1.0 | 2: (0.923423; 5.21508e-02) · 4: (0.998920; 5.74520e-04) · 6: (0.999977; 8.25863e-06) |
| 2.0 | 0.05 | 2: (0.655828; 3.98486e-01) · 4: (0.879209; 6.66369e-02) · 6: (0.965292; 1.31642e-02) |
| 2.0 | 0.1 | 2: (0.659501; 3.90698e-01) · 4: (0.885521; 6.26187e-02) · 6: (0.968652; 1.18702e-02) |
| 2.0 | 0.25 | 2: (0.683646; 3.43942e-01) · 4: (0.918980; 4.27521e-02) · 6: (0.983323; 6.35384e-03) |
| 2.0 | 0.5 | 2: (0.751815; 2.40857e-01) · 4: (0.969677; 1.64883e-02) · 6: (0.996673; 1.39374e-03) |
| 2.0 | 1.0 | 2: (0.883193; 1.03615e-01) · 4: (0.996311; 2.25676e-03) · 6: (0.999858; 6.27382e-05) |

An extrapolated seed is unnecessary for points with an actual certificate; no negative conclusion
about inaccessible sizes follows from this scan.

## Rigor and validation

[METHOD.md](METHOD.md) gives the complete implemented enclosure theorem. In brief, interval
assembly bounds ||A−H||, BigFloat upward accumulation plus explicit dot-product rounding
bounds encloses R=HQ−QΛ and E=Q*Q−I. The polar factor of Q and Weyl's min-max inequality
then enclose every eigenvalue with multiplicity. Interval sums enclose all three log Z values;
outward rational defect bounds feed the unmodified `certified_gap_bound`. Source hashes, Julia
and package versions, dimensions, radii, all eigenvalues and rational seeds are stored in `results/*.json`.

For an exact Hermitian block A, a ≥ ||A−H||₂, r ≥ ||HQ−QΛ||₂ and e ≥ ||Q*Q−I||₂ < 1/2,
the polar factor U obeys ||Q−U||₂ ≤ e. Hence Weyl gives, for every ordered eigenvalue,
**|λᵢ(A)−Λᵢ| ≤ a + r + (||H||F + max|Λᵢ|)e**. This is the theorem implemented;
its norms include the explicit outward rounding allowances proved in METHOD.md.

Certified: RP; complete dimensions; eigenvalue enclosures; partition/seed enclosures; resulting
uniform gap for the stated even-length range. Numerical: seed ranking, ED reference gaps,
ratios, δ(n), purity diagnostics and free-fermion cross-check values. No continuum parameter
region, h→0 limit, odd-size gap, or thermodynamic reference-gap extrapolation is claimed.

The independent audit passed 42 checks: full computational-basis ED versus momentum ED,
the h=0,g=2 `tfim_logZ` formula at six lengths and five temperatures, and exact classical
spectra inside the implemented enclosures. Maximum free-fermion log Z discrepancy: 6.821210263296962e-13.
The fast checker repeats RP, dimension counts, partition bounds and the recurrence. The deep
checker regenerates rigorous eigensystem enclosures and verifies containment in the stored intervals using exact rationals. The stored radii have a factor-two replay margin.
Test commands and measured outcomes are in [README.md](README.md).

## Environment and runtime

Julia 1.12.4; dependencies instantiated and precompiled from the existing Julia project without
modifying its Project or Manifest. OPENBLAS_NUM_THREADS=1 and BLAS.set_num_threads(1).
The approved installer endpoint returned HTTP 403; the official Linux x86_64 archive was used
in /tmp and its SHA-256 checked against the official checksum:
`c57baf178fe140926acb1a25396d482f325af9d7908d9b066d2fbc0d6639985d`.
Sum of per-point measured runtimes: 268.320931 seconds;
audit: 17.481193 seconds. Startup compilation is additional.
All deliverables are confined to this certificate directory. No merge or publication performed.

Final verification: **72/72 fast checks passed in 8.73 s including compilation**;
**30/30 independent spectral replays passed in 248.67 s**; all **1,000 source-hash
comparisons passed**. See [results/validation.json](results/validation.json) for the
machine-readable record and individual certificate fingerprints.

## Independent final audit

A separate reviewer **accepted** the certification argument with **no unresolved theorem
gap**: [INDEPENDENT_REVIEW.md](INDEPENDENT_REVIEW.md). The reviewer confirmed that no
spin-flip parity reduction is used, all translation sectors are included, and the explicit
rounding allowances support the polar/Weyl enclosure theorem.

Additional independent checks passed:

- **1,592 checks** over all 30 stored certificates: character-formula dimensions,
  exact Pauli trace moments, 256-bit interval partition ratios with explicit endpoint
  directions, exact parameter coverage, rational seeds, recurrence and minimum bound.
  See [results/final_audit.json](results/final_audit.json).
- **360 exact rational residual/orthogonality inequalities** on 180 L=6 blocks, covering
  every parameter point. See [results/independent_residual_audit.json](results/independent_residual_audit.json).

The smallest certified bound is exactly **2749898627347949/36028797018963968**,
approximately **0.07632501928667015**, at **g=1, h=1/20**, with **n₀=6, β₀=13/2**.
It applies to every even L≥6 cells=sites. No larger scan was required for this audit.

## Draft PR integration

Upstream subsequently renamed the shared package to `CertifiedGaps`. Certificate imports
now use that package. Exact namespace-only equivalence passed for 29 mathematical/core
sources; original generation sources and hashes are preserved in
[results/original_generation_sources.json](results/original_generation_sources.json), with migration evidence
in [results/namespace_migration.json](results/namespace_migration.json). Original per-point
JSONs and their bounds are unchanged. Complete replay and independent audits run against
the current core. The original independent-review fingerprints refer to the archived
pre-rename sources. This PR changes only the explicitly requested certificate directory
relative to upstream; no core modifications, CI additions, or PR merge are included.

Current-core replay completed successfully before draft PR publication: 72 fast tests
passed in 8.84 s including compilation, all 30 deep spectral replays passed in 251.85 s,
and 1,000 historical/current source-hash comparisons passed. The independent 1,592-check
artifact audit, 360 exact-rational checks and 42 cross-checks also passed on CertifiedGaps.
