# Independent mathematical and code review

Review date: 2026-10-11. This review was performed by a separate reviewer agent,
without modifying the certificate implementation or its existing results.

**Decision: accept the implemented finite-spectrum certification argument.**
No unresolved mathematical gap requiring downgrading these certificates was found.
This decision uses the repository's two-purity theorem and the supported IEEE arithmetic
model stated below. The parent audit separately checks all 30 stored parameter points,
partition bounds, dimensions, hashes, and final bounds. A source review and exact residual
cross-check are not a substitute for that complete stored-artifact replay.

## Hamiltonian, reflection positivity, and symmetry

The integer diagonal is `-L+2*domain_walls-h*(L-2*popcount)` and every single-bit
flip has coefficient `-g`. These are Pauli matrices, not spin operators divided by two.
At L=2 the two directed periodic bonds contribute twice, consistently with the stated
periodic sum. The exact local identity is `-h_bond=Z⊗Z` with coefficient one, and
`-gX-hZ` is real symmetric. `rp_check` checks rational matrix entries exactly.

There is no spin-flip parity reduction. Only translation orbits and all integer
momenta `k=0,...,L-1` are used. Global spin-flip does not commute with the h>0
Hamiltonian; omitting such a sector would have invalidated the result, but the code
does not do so. Reflection exchanges k and -k, and is not used for a further reduction.

An orbit of period p allows precisely the p momenta satisfying kp=0 mod L. Its
Fourier vectors have normalization 1/sqrt(p). Summing over all orbits and momenta
therefore gives exactly 2^L orthonormal vectors. The code asserts this count and
retains every block's full dense spectrum. For a flip of representative b yielding
T^t a, projecting the translated orbit sum gives the implemented coefficient
`-g*sqrt(p_b/p_a)*exp(2πikt/L)`. No numerical approximate basis is assumed unitary.

Interval assembly encloses the exact rational parameters, period ratios, square
roots, and trigonometric phases. Using only the upper triangle to define the binary
Hermitian matrix H is legitimate because the exact block A is Hermitian. The sum of
real and imaginary error radii bounds each entry's modulus error, and the maximum
row sum bounds the spectral norm for the Hermitian difference A-H. Copying that
radius into both symmetric positions is conservative.

## Arithmetic and residual theorem

The analysis does not assume BigFloat by itself certifies a floating calculation.
It includes two distinct errors: interval assembly error, and floating residual-product
error. All binary64 inputs to the BigFloat norm calculation are exact binary rationals;
squares, sums, and square roots are evaluated with RoundUp.

For a d-term complex dot product, expansion into real products and sums, the standard
absolute dot-product error inequality, and
`||abs(B)abs(C)||_F <= ||B||_F ||C||_F` give an error bounded by a small constant
multiple of gamma times the product of Frobenius norms. Here
`gamma=(16d+32)u/(1-(16d+32)u)`, u=2^-53. Complex component comparisons cost at most
a factor two; the factor eight in both allowances safely also covers the final
subtraction and the diagonal multiplication by Lambda. A complex dot product uses
at most 4d real products and 4d real additions, comfortably below the chosen count;
FMA reduces rounding operations. Binary64 matrix products need not use a particular
summation order for this conservative count bound.

Specifically, the code adds

```
8 gamma (||Q||_F² + sqrt(d))
8 gamma ||Q||_F (||H||_F + max|Lambda|)
```

to the computed defect and residual Frobenius norms. These dominate the exact errors
in `Q'*Q-I` and `H*Q-Q*Lambda`, respectively. There is no reversed denominator-rounding
issue: under d<10^6, `(16d+32)*2^-53` and `1-(16d+32)*2^-53` are exact in 256-bit
BigFloat. The subsequent division is rounded upward. The norm and eta additions and
multiplications involve nonnegative upper bounds and are also rounded upward.

The explicit dimension/magnitude restrictions prevent intermediate overflow.
With gradual underflow, even a deliberately loose O(d^3) count of minimum-normal
absolute errors is smaller than the added 2^-900 allowance. The runtime asserts
that flush-to-zero is disabled. This reasoning assumes ordinary IEEE binary64
OpenBLAS operations without fast-math; it is not a proof for arbitrary nonstandard
BLAS implementations or altered floating-point modes. No such alteration appears
in the checked execution setup.

For e >= ||Q'*Q-I||_2 < 1/2, the polar factor U is unitary and
`||Q-U||_2 <= e`, since every singular value satisfies
`|sigma-1|=|sigma²-1|/(sigma+1) <= e`.
The identity

```
HU-U Lambda = (HQ-Q Lambda) + H(U-Q) - (U-Q)Lambda
```

then proves the stated residual bound. Applying the Hermitian min-max/Weyl inequality
to H and U Lambda U* and adding ||A-H|| gives ordered eigenvalue enclosures with
multiplicity. Neither separated eigenvalues nor eigenvector condition estimates are
needed. The stored radius is twice the proved radius; replay containment is checked
with exact rational arithmetic against the fresh proved half-radius.

### Additional independent exact-arithmetic check

I regenerated all six L=6 momentum blocks at every one of the 30 parameter points.
For each dense eigensystem, I converted H, Q, and Lambda to exact complex rational
matrices and computed `HQ-Q Lambda` and `Q'*Q-I` entirely in Rational{BigInt}.
The exact sums of squared moduli were each compared to the square of the implemented
certified upper bound. **All 180 residual checks and 180 orthogonality checks passed**
in 5.857310 seconds, including compilation inside the timed loop. The largest ratios
of exact Frobenius norm to implemented upper bound were 0.00516938495 for residuals
and 0.01703863978 for orthogonality. This is an independent check of the allowances
on actual nontrivial eigensystems, not the mathematical basis of their validity at L=12.

The executable body used for this check is reproduced here; save it to a temporary
Julia file and run from the repository root with the same project/depot as the README:

```julia
include(joinpath(pwd(), "models/tfim_longitudinal/certify.jl"))
using .LongitudinalCertificate, LinearAlgebra
const C = LongitudinalCertificate
q(x) = Rational{BigInt}(x)
cq(z) = complex(q(real(z)), q(imag(z)))
for g in C.GS, h in C.HS, k in 0:5
    H, a = C.block(C.orbits(6), k, g, h; rigorous=true)
    F = eigen(Hermitian(H))
    Q, A, lam = cq.(F.vectors), cq.(H), q.(F.values)
    R, E = A*Q-Q.*transpose(lam), Q'*Q-I
    enc = C.enclosed_eigen(H, a)
    @assert sum(abs2, R) <= q(enc.residual)^2
    @assert sum(abs2, E) <= q(enc.orthogonality)^2
end
```

## Partition ratios, recurrence, and scope

Every spectrum includes multiplicities explicitly. Subtracting the exact common
shift is only a numerically stable evaluation of the same partition function;
the shift need not bound the ground energy. Interval arithmetic encloses both ends
of each log Z. For spatial purity the conservative defect uses
`1-exp(lower(log Z_2n)-2*upper(log Z_n))`; thermal purity uses
`1-exp(lower(log Z_2n(2β))-2*upper(log Z_2n(β)))`.
The implemented interval expressions have these directions, and upper endpoints
are converted exactly to Rational{BigInt}. They are not inward decimal roundings.

Those rational seeds and beta enter the unchanged certified recurrence, whose
admissibility and tail correction are included. The equal-defect 0.12731 condition
is sufficient, not necessary; a successful asymmetric recurrence cannot be rejected
merely because one initial defect exceeds that threshold. The uniform length claim
is the core theorem's every even L>=n0, not just the doubling subsequence. Here one
cell is one site and all chosen n0 are six.

Uniqueness is independently justified for all these finite rings by Perron-Frobenius:
g>0 makes the off-diagonal spin-flip graph irreducible, and a sufficiently large
scalar shift of -H is nonnegative. This argument does not by itself establish a
uniform gap; the certified two-purity theorem supplies that conclusion. No result
is proved for unsampled points, odd rings, or a limit h→0. ED reference gaps and
parameter search rankings remain numerical.

The fast checker validates stored spectral inputs but does not independently prove
them without deep replay. The complete deep replay is therefore essential evidence;
it must not be replaced by the fast test when accepting newly generated certificates.

## Reviewed source fingerprints

Julia runtime independently checked: 1.12.4. SHA-256 at review time:

| Source | SHA-256 |
|---|---|
| certify.jl | `9a45f75b8fba2f90e6b425b66df6b95e892780909ad37e2c8009db0f60d5d9d8` |
| check.jl | `d7920503731dd9422717562f6c06fd0785f65f5fd9ac16876fe06a7e0556ad51` |
| METHOD.md | `cebc4f0903b30448a4a5d1e427f78184cdb92423148ec0e7f10e192153fb4752` |
| ../../julia/src/recurrence/certified.jl | `3ea9426f0a87aed65cbd66f65d442c88437be11056aeb55b386c3b80eb9527e8` |
| ../../outputs/thermal-mpo-gap-scope.tex | `4e77e5d839d79efc4472b6b75a505726b54fb79de5fd4604c27d8927b8d809e7` |

The exact-residual audit is also preserved as
[independent_residual_audit.jl](independent_residual_audit.jl), with its rerun result in
[results/independent_residual_audit.json](results/independent_residual_audit.json):
360/360 exact inequalities passed in 5.911134 seconds. The JSON records the maximum
squared ratios as exact rationals. This is a validation artifact, not a standalone
uniform-gap certificate. Run it with `julia --project=julia
models/tfim_longitudinal/independent_residual_audit.jl`.

The complementary stored-artifact audit reported by the parent passed 1,592 checks
covering all 30 points, using independent character/Ramanujan dimension formulas,
Pauli trace and trace(H²) identities, and 256-bit interval partition ratios:
[results/final_audit.json](results/final_audit.json). These augment the existing full
spectral replay; they do not replace the enclosure theorem reviewed above.

The fingerprint table above refers to the sources as reviewed, before the core package was
renamed from `ThermalMPOGap` to `CertifiedGaps` and the directory moved to
`models/tfim_longitudinal/`; the algorithms are unchanged, and the deep replay
(`DEEP_CHECK=1 check.jl`) re-derives every enclosure with the current core.
