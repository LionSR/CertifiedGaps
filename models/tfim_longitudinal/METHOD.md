# Method: eigenvalue enclosure theorem and arithmetic

## Polar-factor and Weyl eigenvalue enclosure

Let $A$ be an exact Hermitian momentum block and $H$ its `Float64` Hermitian approximation. Interval assembly yields an upper bound $`a \ge \|A - H\|_\infty \ge \|A - H\|_2`$. For the approximate eigenpairs $(Q, \Lambda)$ returned by LAPACK, we treat every entry of $Q$ and $\Lambda$ as an exact binary rational without assuming accuracy from the solver. Let $R = HQ - Q\Lambda$ and $E = Q^\dagger Q - I$, and suppose we have certified bounds $`r \ge \|R\|_2`$ and $`e \ge \|E\|_2`$ with $e < 1/2$. Then the unitary polar factor $U = Q(Q^\dagger Q)^{-1/2}$ satisfies

```math
\|Q - U\|_2 = \max_i |\sigma_i(Q) - 1| \le e,
```

which implies

```math
\|H - U\Lambda U^\dagger\|_2 \le r + (\|H\|_2 + \|\Lambda\|_2)\, e.
```

Applying Weyl's Hermitian eigenvalue perturbation inequality to the sorted eigenvalues gives

```math
|\lambda_i(A) - \Lambda_i| \le \eta := a + r + (\|H\|_F + \max_j |\Lambda_j|)\, e
```

for every ordered eigenvalue with multiplicity, without requiring spectral gap separation between adjacent eigenvalues.

To ensure platform-independent replay across different LAPACK implementations, the certificate stores the doubled radius **$2\eta$** (multiplication by 2 is exact in `BigFloat`), and `check.jl --full` verifies in exact rational arithmetic that every freshly computed interval $[\Lambda_i^{\text{new}} - \eta^{\text{new}}, \Lambda_i^{\text{new}} + \eta^{\text{new}}]$ is contained inside the stored interval $[\Lambda_i^{\text{stored}} - 2\eta^{\text{stored}}, \Lambda_i^{\text{stored}} + 2\eta^{\text{stored}}]$.

## Implemented binary64 rounding allowances

The function `enclosed_eigen` computes the residual and orthogonality defect matrices in `Float64`, and then accumulates all squared absolute entries of $R$, $E$, $H$, and $Q$ in 256-bit `BigFloat` with `RoundUp`. To account for the floating-point rounding error committed during the `Float64` matrix multiplications $Q^\dagger Q - I$ and $HQ - Q\Lambda$, let $d$ be the block dimension, $u = 2^{-53}$, $m = 16d + 32$, and $\gamma_m = mu / (1 - mu)$. The implemented bounds are

```math
e = \|\mathrm{fl}(Q^\dagger Q - I)\|_F + 8\gamma_m (\|Q\|_F^2 + \sqrt{d}) + 2^{-900},
```
```math
r = \|\mathrm{fl}(HQ - Q\Lambda)\|_F + 8\gamma_m \|Q\|_F (\|H\|_F + \max_j |\Lambda_j|) + 2^{-900}.
```

These bounds follow from the standard IEEE binary64 dot-product bound $|\mathrm{fl}(\sum a_j b_j) - \sum a_j b_j| \le \gamma_m \sum |a_j b_j|$, expansion of complex products into real operations, and the submultiplicativity $`\||B|\,|C|\|_F \le \|B\|_F \|C\|_F`$. The code asserts $d < 10^6$ and entry magnitudes below $10^6$ so that no overflow can occur, and the additive $2^{-900}$ term bounds gradual-underflow errors across $O(d^3)$ operations.

## Exact translation-sector reduction and interval assembly

Translation $T$ acts by cyclic rotation of the $L$-bit spin configuration. For each bitstring orbit representative $s$ of minimal period $p$ and momentum index $`k \in \{0, \dots, L-1\}`$ satisfying $kp \equiv 0 \pmod L$, the normalized translation state is

```math
|s; k\rangle = \frac{1}{\sqrt{p}} \sum_{r=0}^{p-1} e^{-2\pi i k r / L} |T^r s\rangle.
```

An orbit of period $p$ contributes $p$ states across all momenta $k$, and the code verifies that the sector dimensions sum to $\sum p = 2^L$. If flipping a single spin in representative $b$ yields $T^t a$, its off-diagonal matrix element is $`-g\sqrt{p_b / p_a}\, e^{2\pi i k t / L}`$; the diagonal entry is $-L + 2(\text{domain walls}) - h(L - 2\mathrm{popcount}(s))$ (with both periodic bonds included at $L = 2$).

All parameters, orbit periods, and diagonal counts are exact integers or rationals. `IntervalArithmetic` encloses every square root, $\sin(\pi x)$, and $\cos(\pi x)$ entry, converts the midpoint to a Hermitian `Float64` matrix $H$, and bounds $`\|A - H\|_\infty`$ via the maximum row sum of component radii. Every momentum sector $k = 0, \dots, L-1$ is diagonalized independently without assuming reflection or spin-flip multiplicities.

## Partition functions and rational seed defects

All stored eigenvalues and enclosure radii are serialized as exact rational strings. With $c = \min_i \Lambda_i$ chosen to stabilize the exponentials, `IntervalArithmetic` encloses

```math
\log Z_L(\beta) = -\beta c + \log \sum_i e^{-\beta(\lambda_i - c)}
```

across all enclosed eigenvalues $\lambda_i \in [\Lambda_i - 2\eta_i, \Lambda_i + 2\eta_i]$. From the enclosed intervals for $\log Z_n(\beta_0)$, $\log Z_{2n}(\beta_0)$, and $\log Z_{2n}(2\beta_0)$, the upper endpoints of $1 - P_{\mathrm{s}}(n, \beta_0)$ and $1 - P_\tau(2n, \beta_0)$ are converted **exactly** to `Rational{BigInt}` strings and passed to `CertifiedGaps.certified_gap_bound`.

The fast checker repeats the RP check, the dimension counts, the interval partition sums, and the
recurrence. `check.jl --full` additionally regenerates the dense eigenvalue enclosures and verifies
their containment in the stored spectral inputs; no eigenvectors are trusted from a JSON file.
