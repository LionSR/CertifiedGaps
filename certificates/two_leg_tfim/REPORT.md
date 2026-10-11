# Two-leg TFIM: verification method and results

**Result: g=3 is the smallest certified field on the tested grid.** The seed n0=4 rungs, β0=3/2 proves Δ_L ≥ 0.636386863 for every even L≥4 rungs (N=2L spin sites, N divisible by four and N≥8). The exact returned binary bound is recorded below. Fields 2 and 5/2 fail the numerical grid search; g=4 passes numerically but was not rerun through verification. Only the best seed at the smallest successful field was verified.


## Model and scope

One cell is a rung, with two spin-1/2 sites and dimension four. In Pauli normalization,
H_L = −Σ_j(Z_{j,a}Z_{j+1,a}+Z_{j,b}Z_{j+1,b}+Z_{j,a}Z_{j,b})
−gΣ_j(X_{j,a}+X_{j,b}), with periodic rung index. The length-two convention includes
both periodic directed bonds, as in the core package. Fields are exactly 2, 5/2, 3, 4.
For every g>0, the σᶻ-basis off-diagonal graph is the connected hypercube with negative
edge weights −g. Perron–Frobenius applied to cI−H proves a unique finite-ring ground
state. This alone does not establish a uniform thermodynamic gap.

The remembered critical field 1.6 is inaccurate for this normalization. Xavier, Pereira,
Nunes and Plascak, [arXiv:2111.11404, Eq. (1), Table I](https://arxiv.org/html/2111.11404),
set J_perp^x=0, J_perp^z=1; at J4=0 their DMRG extrapolation is g_c=1.83214(5).
This numerical literature result places all requested fields on the paramagnetic side;
it is not used in the certificate.

## Exact RP and symmetries

In the real σᶻ cell basis (↑↑,↑↓,↓↑,↓↓), set
Z=diag(1,−1), X=[[0,1],[1,0]], C1=Z⊗I, C2=I⊗Z. Then
−h_bond=C1⊗C1+C2⊗C2, with exact nonnegative rational coefficients c1=c2=1.
The onsite matrix is o=−Z⊗Z−g(X⊗I+I⊗X), which is real symmetric.
`rp_identity` checks the four-spin matrix identity with Rational{BigInt}, independently
constructing the bond using four 2×2 Kronecker factors. No floating-point RP test is
substituted for this identity. As a single rung cell, reflection exchanges identical
left and right leg bonds and preserves the onsite term; there is no nonpalindromic
blocking issue.

For ED only, a real Hadamard rotation on every spin diagonalizes the field. The parity
of the number of down spins is conserved; pair flips realize the ZZ couplings. Leg
exchange S, rung reflection R:j↦−j, and (for even L) half translation T:j↦j+L/2
are commuting involutions preserving H. For each parity and real character χ, project
each permutation orbit: |O,χ⟩=Σ_{s∈O}χ(s)|s⟩/√|O| if its stabilizer has trivial χ,
and omit it otherwise. Distinct orbits have disjoint support; character orthogonality
resolves each orbit. Thus these are exact orthonormal invariant blocks, exhaustively
counted by the asserted sum of block dimensions 4^L.

Sector records are ordered by parity 0,1 and integer character labels. Character bits select leg exchange, half translation, and reflection (the half-translation bit is omitted for odd L).

All block spectra are obtained by dense real symmetric LAPACK diagonalization, without
truncation or missing multiplicities. Exact block entries are A_ij/2√(m_j/m_i),
where A has integer entries and m_i are orbit sizes. Exact symmetry is checked by
A_ij m_j=A_ji m_i. The small-ring tests compare with an independently assembled full
Hamiltonian. The theorem in the scope note explicitly requires even n0. The search
includes n0=2,3,4, but only 2 and 4 can be selected for certification. These are all
even seeds whose doubled rings have at most ten rungs. The ten-rung odd-seed diagnostic
n0=5 is not computed and no claim of exhaustive odd-seed search is made.

## Implemented enclosure theorem

Let B be the exact real symmetric block, H its binary floating-point approximation,
Q the square matrix of computed eigenvectors, and Λ the diagonal matrix of computed
real eigenvalues. The verifier establishes
||B−H||₂≤d, ||HQ−QΛ||₂≤r, ||QᵀQ−I||₂≤e<1.
For the polar factor U=Q(QᵀQ)^(−1/2), U is orthogonal and
||U−Q||₂≤a=e/(1+√(1−e)). Hence
||UᵀHU−Λ||₂≤r+(||H||₂+||Λ||₂)a.
Weyl's ordered-eigenvalue inequality and the triangle inequality give, for every
ordered eigenvalue with multiplicity,
|λ_i(B)−Λ_i|≤η=d+r+(||H||_F+||Λ||₂)a.
This proves a whole-spectrum enclosure, unlike a residual-only bound on individual
Ritz values which would not establish multiplicity coverage.

Implementation details: each algebraic B entry is enclosed with directed 256-bit MPFR
square root/division/multiplication. The maximum row sum of construction errors bounds
d because B−H is symmetric. The stored floating-point R and E are recomputed from the
returned full eigenvector matrix. The standard floating-point dot-product error bound
|fl(AB)−AB|≤γ_k |A||B|, with u=2^−53 and γ_k=ku/(1−ku), is applied with the
conservative k=2m+8 including multiplication and subtraction. Thus the implementation
adds γ_k(||H||_F||Q||_F+||Q||_F||Λ||₂) to the residual bound and
γ_k(||Q||_F²+√m) to the orthogonality bound. Stored matrix norms are bounded by
m times their maximum absolute entry. The scalar Frobenius norms are accumulated with
BigFloat rounding upward, as are all combined upper bounds. An additive
16m³·2^−1074 covers gradual-underflow errors. This assumes IEEE binary64 arithmetic
with gradual underflow and conventional BLAS dot products, including possible FMA;
no assumption of exact floating-point eigenvalues is made.

Directed MPFR exponential, summation and logarithm enclose every log Z from the full
spectrum with λ_i∈[Λ_i−η,Λ_i+η]. In particular lower bounds are required for the
numerators, not merely upper bounds on all three log Z values. The lower logarithmic
purities are l_s=lo(log Z_2n(β))−2hi(log Z_n(β)) and
l_t=lo(log Z_2n(2β))−2hi(log Z_2n(β)). Upward-rounded −expm1(l_s), −expm1(l_t)
are converted exactly from their binary BigFloat values to rational strings. These
rationals and the exact dyadic β are passed to `certified_gap_bound` unmodified.
The returned `.bound` is the uniform lower bound, not the finite-iteration upper rate.

## Results

All scan values below are **numerical, not certified**. The grid is β=1/8,2/8,…,32; n0=2,3,4.

| g | best eligible n0 | β | numerical bound | ED gap (8 rungs) | spectral runtime (s) |
|---|---:|---:|---:|---:|---:|
| 2//1 | — | — | no seed | 0.528639460 | 89.62 |
| 5//2 | — | — | no seed | 1.579207409 | 134.27 |
| 3//1 | 4 | 1.500 | 0.636872819 | 2.661781588 | 90.40 |
| 4//1 | 4 | 1.125 | 2.142199085 | 4.763362747 | 136.61 |

Spectral runtimes sum the measured dense-spectrum calls after orbit construction; they exclude Julia startup, the β-grid scan, and JSON I/O. Verification runtime includes spectra, residual validation, and seed construction, excluding startup and final JSON I/O. Each completed run is below the 45-minute limit.


### Purity and ground-energy finite-size diagnostics

| g | n | max_β min(Ps,Pτ) | δ(n)=E0(2n)−2E0(n) |
|---|---:|---:|---:|
| 2//1 | 2 | 0.454349841 | 0.753150185 |
| 2//1 | 3 | 0.526080292 | 0.393739787 |
| 2//1 | 4 | 0.564897524 | 0.262460357 |
| 5//2 | 2 | 0.532084570 | 0.717774565 |
| 5//2 | 3 | 0.683199040 | 0.292564344 |
| 5//2 | 4 | 0.780193496 | 0.142802875 |
| 3//1 | 2 | 0.602242471 | 0.634477163 |
| 3//1 | 3 | 0.797921242 | 0.200727919 |
| 3//1 | 4 | 0.893759824 | 0.073955799 |
| 4//1 | 2 | 0.722816094 | 0.492692762 |
| 4//1 | 3 | 0.905693809 | 0.106397922 |
| 4//1 | 4 | 0.969496099 | 0.026210081 |

For the unsuccessful fields, a two-point fit δ(n)=δ(4)exp(−(n−4)/ξδ) to n=3,4, combined with an extensive dilute-excitation extrapolation of the eight-rung partition function, predicts the following first even seeds. These are **uncontrolled estimates**, not computed larger-ring spectra or certificates; near-critical finite-size effects can bias them substantially. `extrapolate.jl` records the explicit ansatz and grid.

- g=2//1: ξδ=2.465543 rungs; predicted (n0,β)=(10,8.250), requiring 20 rungs for direct verification.
- g=5//2: ξδ=1.394273 rungs; predicted (n0,β)=(6,2.500), requiring 12 rungs for direct verification.

### Dense sector dimensions

- L=2: 5, 1, 1, 1, 2, 2, 2, 2; sum=16=4^L. Zero-dimensional sectors are omitted.
- L=3: 14, 6, 6, 6, 10, 10, 6, 6; sum=64=4^L. Zero-dimensional sectors are omitted.
- L=4: 34, 18, 14, 14, 10, 10, 14, 14, 24, 24, 16, 16, 8, 8, 16, 16; sum=256=4^L. Zero-dimensional sectors are omitted.
- L=6: 320, 256, 256, 256, 240, 240, 240, 240, 272, 272, 272, 272, 240, 240, 240, 240; sum=4096=4^L. Zero-dimensional sectors are omitted.
- L=8: 4384, 4128, 4064, 4064, 4000, 4000, 4064, 4064, 4224, 4224, 4096, 4096, 3968, 3968, 4096, 4096; sum=65536=4^L. Zero-dimensional sectors are omitted.

### Certified results

- **g=3//1: Δ_L ≥ 0.636386863**, even L >= 4 rungs (2L spin sites, multiples of four, >= 8 sites). Seed (n0,β)=(4,3//2); verification 645.83 s; bound / 8-rung ED gap = 0.239083.
  Returned `.bound` = `0.6363868632722607` (exact binary rational `179126977518527//281474976710656`); the displayed decimal lower bound is rounded downward.
  Exact rational seeds and per-sector eigenvalue radii are in `results/certificate_g3_1.json`.

Certified: the exact RP identity, whole-spectrum eigenvalue enclosures for the selected parameter/seed, outward log-partition and rational defect bounds, and the resulting uniform gap under the repository's two-purity theorem. Numerical: seed optimization, all nonverified scan spectra and finite-size gaps, the literature critical point, and comparisons to finite-ring reference gaps. The ED reference gap is not claimed to be an infinite-volume limit.


## Validation and environment

Julia 1.12.4 was verified; `Pkg.instantiate()` and `Pkg.precompile()` succeeded with the existing `julia/` project. The juliaup endpoint returned HTTP 403, so the official versioned Julia tarball was used. A writable temporary depot resolved the read-only home directory. `results/environment.json` records package versions, BLAS configuration, and hashes of the project, manifest, repository instructions, and scope note.

`julia --project=julia certificates/two_leg_tfim/run_tests.jl` passed **69 tests in 19.93 seconds**: 10 exact-RP/symmetry checks, 57 stored-seed/recurrence/provenance checks, and 2 enclosure checks on an exactly soluble matrix. `results/tests.json` records the measured runtime and source hashes. The fast checker recomputes log Z and rational seeds from stored eigenvalue radii, not the expensive residuals; rerun `certify.jl 3//1` to reproduce those.

The completed verification kernel took 645.83 seconds; including final JSON serialization the script reported 652.22 seconds. Maximum eigenvalue enclosure radii were approximately 1.50e-10 for four rungs and 3.94e-5 for eight rungs. The outward rational defect bounds have decimal approximations εs≈0.11240836635 and εt≈0.05161523886; exact fractions are in the certificate.

No numerical or dependency blocker remains. Search scope is n0≤4, with doubled rings through eight rungs. This covers every even seed permitted by the scope-note theorem whose doubled ring fits the requested ten-rung cap. The odd n0=5 / ten-rung diagnostic was not computed. No conclusion of impossibility at larger sizes or other temperatures is claimed.
