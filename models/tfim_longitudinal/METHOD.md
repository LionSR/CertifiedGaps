# Enclosure theorem and arithmetic

Let A be an exact Hermitian momentum block and H its stored binary64 Hermitian
approximation. Interval assembly gives a >= ||A-H||_infinity >= ||A-H||_2.
For the complete dense eigensystem returned by LAPACK, treat every entry of Q and
Lambda as an exact binary rational; do not assume LAPACK's output is exact.
Let R=HQ-Q Lambda and E=Q*Q-I. Suppose r >= ||R||_2 and e >= ||E||_2, with e<1/2.
The polar factor U=Q(Q*Q)^(-1/2) is unitary and

    ||Q-U||_2 = max_i |sigma_i(Q)-1| <= e.

Consequently

    ||H-U Lambda U*||_2 <= r + (||H||_2+||Lambda||_2)e.

Weyl's Hermitian eigenvalue perturbation inequality, applied in sorted order, yields

    |lambda_i(A)-lambda_i(Lambda)| <= a+r+(||H||_F+max|Lambda|)e = eta.

For stable independent replay, the stored radius is **2 eta**; multiplication by two
is exact for these finite BigFloat values. Replay checks that newly certified intervals
of radius eta fit inside the stored intervals, not bitwise equality of eigensolves.

This matches every eigenvalue with multiplicity, including degenerate ones. There is
no spectral separation assumption and no reliance on small residuals of individual vectors.
The polar identity and the displayed triangle inequality prove the needed residual theorem;
Weyl follows directly from the min-max principle.

## Implemented rounding allowances

`enclosed_eigen` computes the residual and defect matrices in binary64, then accumulates
**all squared absolute entries** of these matrices, H and Q in 256-bit BigFloat using
RoundUp. Each binary64 input is exactly representable in BigFloat. These sums give
outward Frobenius norms of the computed matrices, not of the unrounded products.
The difference is covered explicitly. With d the block dimension, u=2^-53,
m=16d+32, gamma=m*u/(1-m*u), the implemented allowances are

    e = ||fl(Q*Q-I)||_F + 8 gamma (||Q||_F^2 + sqrt(d)) + 2^-900
    r = ||fl(HQ-Q Lambda)||_F
          + 8 gamma ||Q||_F (||H||_F+max|Lambda|) + 2^-900.

They follow from the standard dot-product rounding bound
|fl(sum a_j b_j)-sum a_j b_j| <= gamma_m sum |a_j b_j|, expansion of
complex multiplication into real products and sums, and
|| |B| |C| ||_F <= ||B||_F ||C||_F. Splitting real and imaginary parts
costs at most a factor two per complex modulus comparison; the factor eight and
m=16d+32 also cover the final subtraction and diagonal scaling. These deliberately
loose allowances hold for arbitrary summation order with at most these operation
counts (including FMA); they are not empirical estimates of BLAS accuracy.
The code asserts d<10^6 and magnitudes of H,Q,Lambda below 10^6. Products cannot
overflow. With gradual underflow, all additive underflow errors, even bounded by
O(d^3) binary64 minimum-normal values, are below the extra 2^-900 allowance.
RoundUp also applies to gamma, e, r and eta. The run uses single-threaded OpenBLAS,
IEEE binary64, no fast-math, and gradual underflow. Fast-math/flush-to-zero modes are
outside the supported arithmetic environment.

## Exact symmetry reduction and interval assembly

Translation T cyclically rotates the bit word. For each orbit representative s of
minimal period p, and k with kp=0 modulo L, the exact vector is

    |s;k> = p^(-1/2) sum_{r=0}^{p-1} exp(-2 pi i k r/L) |T^r s>.

These vectors are orthonormal. All allowed k and all orbits form a complete basis:
an orbit of period p contributes p vectors, so the sum of sector dimensions is
sum p=2^L, asserted for every ring. If a single spin flip of representative b gives
T^t a, its block entry is -g sqrt(p_b/p_a) exp(2 pi i k t/L). The diagonal is
-L+2(number of domain walls)-h(L-2 popcount(s)). Both bonds are included at L=2.

All parameters, orbit indices, periods and diagonal integers are exact. IntervalArithmetic
outward-encloses each rational, square root, sinpi and cospi and accumulates every
entry. The upper triangle is converted to a Hermitian binary64 matrix. The sum of
real and imaginary component error radii in each row bounds ||A-H||_infinity.
No claim rests on a floating-point change of basis being unitary.

Every momentum sector is solved densely; k and -k are both computed. Spin-flip is
not a symmetry for h>0. Reflection exchanges k with -k, so parity is not an independent
label in a generic k block. We do not split the k=0,pi blocks further by reflection.
This sacrifices a constant-factor optimization and avoids any inferred multiplicities.

## Partition functions and rational defects

All stored eigenvalues and radii are exact rational strings. For an exact common shift c,
interval evaluation encloses

    log Z(beta) = -beta c + log sum_i exp[-beta(lambda_i-c)]

using every enclosed eigenvalue. The shift is the smallest computed binary eigenvalue,
not a claimed ground-energy lower bound. It only stabilizes the exponentials.
Compute the three log-partition intervals a,c,d for (n,beta),(2n,beta),(2n,2beta).
Upper endpoints of 1-exp(c-2a) and 1-exp(d-2c) are converted **exactly** from binary64
to rational strings and passed with exact beta to `certified_gap_bound`.
Both lower and upper partition endpoints are required. Finally its `.bound` is recorded
also as an exact binary rational, so the displayed decimal is only a readable rendering.

The fast checker repeats RP, dimension counts, interval partition sums and the recurrence.
`DEEP_CHECK=1` additionally regenerates the dense eigenvalue enclosures, verifying containment in the
stored spectral inputs; no eigenvectors are trusted from a JSON file.
