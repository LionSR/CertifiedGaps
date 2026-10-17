# Method: two-leg Heisenberg ladder certificate

The first part of this document describes the verification pipeline implemented in `src/`. The second part proves that a naive perturbation expansion around independent rungs fails at every seed size and derives the local coupling continuation used in the certificate.

## Certificate pipeline

### Complete verified spectra (`src/spectra.jl`)

The Hilbert space is decomposed by down-spin number, rung momentum, leg exchange, and (at zero magnetization) spin inversion using exact permutation orbits and stabilizers. Opposite magnetizations and conjugate momenta carry exact multiplicities, and the block dimensions are verified to sum to $4^L$. Block entries are enclosed in interval arithmetic and centered to a `Float64` matrix $H$ with $`\|H_{\text{exact}} - H\|_2 \le \text{herror}`$ (bounded via the maximum absolute row sum).

LAPACK provides tentative eigenpairs $(Q, \Lambda)$. Using interval Frobenius bounds $`r \ge \|HQ - Q\Lambda\|_F`$ and $`\eta \ge \|Q^\dagger Q - I\|_F < 1`$—both including the explicit floating-point rounding allowance $`\gamma_{8n+32}`$ (with unit roundoff $u = 2^{-53}$) and a $2^{-1000}$ underflow allowance—the unitary polar factor $U$ satisfies $`\|Q - U\|_2 \le d = \eta / (1 + \sqrt{1 - \eta})`$. Weyl's inequality then gives the simultaneous eigenvalue enclosure radius

```math
\rho = r + (\|H\|_F + \max_i |\Lambda_i|)\, d + \text{herror}.
```

Every eigenvalue ball, weighted by its sector multiplicity, is summed with outward-directed rounding to enclose $\log Z_4(\beta_0)$, $\log Z_8(\beta_0)$, and $\log Z_8(2\beta_0)$.

### Seed classification (`src/suite.jl`)

For each coupling, the suite computes the verified spectra of $2$-, $4$-, and $8$-rung rings and classifies all 128 grid seeds $`n_0 \in \{2, 4\}`$, $\beta_0 = k/4$ ($k = 1, \dots, 64$):

- `certified`: the upper defect bounds $\varepsilon_{\mathrm{s}}, \varepsilon_\tau$ (rounded upward to multiples of $10^{-8}$) close `certified_gap_bound`, and the resulting gap lower bound is rounded downward to a multiple of $10^{-9}$;
- `rejected_by_recurrence`: the lower defect bounds (rounded downward) carry an exact rational divergence witness. With $z = 1 - (1 - \varepsilon_\tau)(1 - \varepsilon_{\mathrm{s}})^2$ and the squaring map $f(z) = z^2 / [2(1 - z)^2]$, once $z \ge 1/2$ (or the intermediate thermal step reaches $1/2$), both components at the next step are at least $1/2$, so the order-preserving recurrence cannot converge to zero (intermediate lower bounds are rounded downward to $10^{-24}$);
- `unresolved`: neither condition holds (this does not occur for any grid seed in the certificate).

### All-$`\beta`$ exclusion ($J_\perp = 2, 3/2$)

Given ground-energy enclosures $a_L \le E_0(L) \le b_L$ from the eigenvalue balls, the shifted log-partition function $`F_L(\beta) = \log \mathrm{Tr}\, e^{-\beta(H_L - a_L)}`$ is monotonically decreasing in $\beta$. On any temperature interval $[l, u]$:

```math
\log P_{\mathrm{s}}(n, \beta) \le F_{2n}(l) - 2F_n(u) - \beta(a_{2n} - 2a_n),
```

evaluated over the interval without assuming monotonicity of $P_{\mathrm{s}}$. Because $P_\tau(2n, \beta)$ is monotonically increasing in $\beta$ ($\frac{d}{d\beta}\log P_\tau = 2(\langle H \rangle_\beta - \langle H \rangle_{2\beta}) \ge 0$), we also have $P_\tau(2n, \beta) \le P_\tau(2n, u)$ on $[l, u]$. Starting from $[0, B]$, intervals are recursively bisected until every sub-interval has an exact recurrence rejection witness. For the low-temperature tail $\beta \ge B$, a positive ground-energy mismatch $\delta = a_{2n} - 2b_n > 0$ yields $\log P_{\mathrm{s}}(n, \beta) \le F_{2n}(B) - B\delta$. Running `check.jl --full` verifies contiguous interval coverage, every window witness, and the tail bound.

### Coupling continuation

For $|J_\perp - J_*| \le \delta$, the centered rung sum has operator norm at most $L/2$, so

```math
\varepsilon_{\mathrm{s}}(\delta) = 1 - (1 - \varepsilon_{\mathrm{s}}) e^{-2\beta_0 n_0 \delta}, \qquad \varepsilon_\tau(\delta) = 1 - (1 - \varepsilon_\tau) e^{-4\beta_0 n_0 \delta}
```

are valid defect upper bounds across the entire interval $|J_\perp - J_*| \le \delta$ (derived below). The suite stores the resulting uniform gap lower bound for $\delta = 10^{-4}$.

## Why the naive strong-rung norm sandwich fails at every length

This section records a rigorous negative result for the naive independent-rung operator-norm perturbation bound: for every even seed length $n \ge 2$, every $\beta > 0$, and every $0 < J_\perp \le 4$, bounding $\log Z_L(\beta)$ around the independent-rung limit via Weyl's inequality cannot produce an admissible seed. Increasing $n$ cannot rescue the bound because the extensive norm error destroys the bulk cancellation between purities.

Let $H_0$ contain only the rung interactions, let $V$ contain the $2L$ leg bonds, and let $z(\beta) = e^{3\beta J_\perp/4} + 3e^{-\beta J_\perp/4}$ be the single-rung partition function. Each spin-1/2 Heisenberg leg bond $h$ has eigenvalues $-3/4$ and $1/4$. Centering $h$ by $+I/4$ gives $`\|h + I/4\|_2 = 1/2`$, so the sum of $2L$ leg bonds satisfies

```math
\|V + \tfrac{L}{2} I\|_2 \le L, \qquad \bigl|\log Z_L(\beta) - \bigl(L \log z(\beta) + \tfrac{1}{2}\beta L\bigr)\bigr| \le \beta L.
```

While this Weyl sandwich holds for all $L$, applying it independently to $Z_n(\beta)$, $Z_{2n}(\beta)$, and $Z_{2n}(2\beta)$ yields the defect upper bounds

```math
\varepsilon_{\mathrm{s}} = 1 - e^{-4\beta n}, \qquad \varepsilon_\tau = 1 - e^{-8\beta n} p(\beta)^{2n}, \qquad p(\beta) = \frac{1 + 3q^2}{(1 + 3q)^2}, \quad q = e^{-\beta J_\perp}.
```

These bounds never enter the admissible region for any $n \ge 2$, $\beta > 0$, and $0 < J_\perp \le 4$:

- If $\beta n \ge 1/8$, then $(1 - \varepsilon_{\mathrm{s}})^2 = e^{-8\beta n} \le e^{-1} < 1/2$.
- If $\beta n < 1/8$, then $\beta J_\perp \le 4\beta < 1/(2n) \le 1/4$. Hence $3/4 < 1 - \beta J_\perp \le q \le 1$, which implies $p(\beta) \le 4 / (1 + 9/4)^2 = 64/169 < 1/2$, and therefore $1 - \varepsilon_\tau < 1/2$.

In both cases, $A = (1 - \varepsilon_{\mathrm{s}})^2(1 - \varepsilon_\tau) < 1/2$, so $z = 1 - A > 1/2$ and $f(z) = z^2 / [2(1 - z)^2] > 1/2$. Consequently, the recurrence diverges above $1/2$ regardless of $n$ or $`\beta`$.

## Local coupling continuation around a certified seed

At a certified coupling $`J_*`$, write $`H_L(J_\perp) - H_L(J_*) = (J_\perp - J_*) W_L`$, where $`W_L`$ is the sum of $`L`$ rung operators. Centering each rung operator by $`+I/4`$ gives $`\|W_L + \tfrac{L}{4} I\|_2 \le L/2`$. For $`|J_\perp - J_*| \le \delta`$, the scalar shift cancels in both purities and Weyl's inequality bounds the variation of each centered log-partition function by $\frac{1}{2}\beta \delta L$. Thus, given certified defects $(\varepsilon_{\mathrm{s}}, \varepsilon_\tau)$ at $`(n_0, \beta_0, J_*)`$, valid upper defect bounds across the entire interval $`|J_\perp - J_*| \le \delta`$ are

```math
\varepsilon_{\mathrm{s}}(\delta) = 1 - (1 - \varepsilon_{\mathrm{s}}) e^{-2\beta_0 n_0 \delta}, \qquad \varepsilon_\tau(\delta) = 1 - (1 - \varepsilon_\tau) e^{-4\beta_0 n_0 \delta}.
```

Evaluating these bounds with outward rounding and running the doubling recurrence certifies the full interval $[J_* - \delta, J_* + \delta]$ without assuming monotonicity in $J_\perp$. Both `check.jl` modes verify this continuation for $\delta = 10^{-4}$.

## Smaller-seed identity

Directly from the definition of the purities,

```math
P_\tau(2n, \beta) = \frac{P_{\mathrm{s}}(n, 2\beta)\, P_\tau(n, \beta)^2}{P_{\mathrm{s}}(n, \beta)^2}.
```

Because reflection positivity implies $P_{\mathrm{s}}(n, \beta) \le 1$ for even $n$, the numerator $P_{\mathrm{s}}(n, 2\beta) P_\tau(n, \beta)^2$ is a valid lower bound on $P_\tau(2n, \beta)$. However, evaluating $P_{\mathrm{s}}(n, 2\beta)$ still requires $Z_{2n}(2\beta)$ unless an independent lower bound on $P_{\mathrm{s}}(n, 2\beta)$ is available. We therefore compute $Z_{2n}(2\beta)$ directly from the complete 8-rung spectrum.
