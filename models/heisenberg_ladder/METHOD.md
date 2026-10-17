# Method: two-leg Heisenberg ladder certificate

The first part of this document describes the frame and the verification pipeline (the shared core plus `src/`). The second part proves that a naive perturbation expansion around independent rungs fails at every seed size and derives the local coupling continuation used in the certificate.

## Certificate pipeline

### Frame and reflection positivity (`src/model.jl`)

The ladder is a core `RingModel` with one rung ($d = 4$) per cell in the leg-rotated frame $U_L = \bigotimes_j (1 \otimes u)$, $u = \exp(i\pi S^y)$: the leg bonds $\mathbf{S} \cdot \mathbf{S}$ are unchanged (both endpoints are rotated), the rung term becomes $o = -J_\perp(X \otimes X + A \otimes A + Z \otimes Z)$, and $S^z_{\text{tot}}$ becomes the additive cell charge $z_1 - z_2$. With the cell reflection $F = \sigma^y \otimes \sigma^y$, which commutes with $o$, the bond satisfies $-h = \sum_\alpha C_\alpha \otimes F C_\alpha F$ over the six real single-spin factors, the reflected RP form of the scope note. Unlike the earlier checkerboard frame (alternating legs, conserved charge only up to a staggering), this frame is translation invariant with a uniform U(1) charge, so the generic core blocks apply. `check_frame` verifies the frame identities exactly.

### Complete verified spectra (core `VerifiedSpectra`)

`verified_ring_spectrum` splits each ring exactly by charge, momentum, the parity of the on-site symmetry leg exchange $`U\,\mathrm{SWAP}\,U^\mathsf{T}`$ and, in the self-conjugate sector $S^z = 0$, the parity of the flip $F^{\otimes L}$; sectors $\pm S^z$ and momenta $\pm k$ are built once with exact multiplicities, and the block dimensions sum to $4^L$ (91 blocks, largest 715, at 8 rungs). Block entries are enclosed in interval arithmetic; the exactly Hermitian Float64 center and an exact Frobenius bound on its distance to the true block go to `verified_eigen`, which bounds the residual and orthogonality of the LAPACK eigenvectors in Float64 with the a priori allowances $\gamma_{N+8}$ (real blocks) or $\gamma_{2N+8}$ (complex blocks, four real products) on Frobenius norms, takes the smaller of the Ostrowski and polar-factor radii, and adds the block error by Weyl's inequality. The largest radius at 8 rungs is $8.6 \times 10^{-9}$ (the earlier local enclosure gave $4.4 \times 10^{-8}$). `log_partition_enclosure` sums every eigenvalue ball with its multiplicity in outward-rounded interval arithmetic.

### Seed classification (`src/suite.jl`)

For each coupling, `seed_bounds` (reusing the three spectra) encloses both defects of all 128 grid seeds $`n_0 \in \{2, 4\}`$, $\beta_0 = k/4$ ($k = 1, \dots, 64$), with upper ends rounded up to multiples of $10^{-15}$:

- `certified`: the upper defect bounds close `certified_gap_bound`; the best such seed is passed to `certify_gap`, which also runs the exact rational chain;
- `rejected_by_recurrence`: the lower defect bounds (rounded downward to multiples of $10^{-15}$) carry an exact rational divergence witness. With $z = 1 - (1 - \varepsilon_\tau)(1 - \varepsilon_{\mathrm{s}})^2$ and the squaring map $f(z) = z^2 / [2(1 - z)^2]$, once $z \ge 1/2$ (or the intermediate thermal step reaches $1/2$), both components at the next step are at least $1/2$, so the order-preserving recurrence cannot converge to zero (intermediate lower bounds are rounded downward to $10^{-24}$);
- `unresolved`: neither condition holds (this does not occur for any grid seed in the certificate).

### All-$\beta$ exclusion ($J_\perp = 2, 3/2$)

Given ground-energy enclosures $a_L \le E_0(L) \le b_L$ from the eigenvalue balls, the shifted log-partition function $F_L(\beta) = \log \mathrm{Tr} e^{-\beta(H_L - a_L)}$ is monotonically decreasing in $\beta$. On any temperature interval $[l, u]$:

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

In both cases, $A = (1 - \varepsilon_{\mathrm{s}})^2(1 - \varepsilon_\tau) < 1/2$, so $z = 1 - A > 1/2$ and $f(z) = z^2 / [2(1 - z)^2] > 1/2$. Consequently, the recurrence diverges above $1/2$ regardless of $n$ or $\beta$.

## Local coupling continuation around a certified seed

At a certified coupling $J_*$, write $H_L(J_\perp) - H_L(J_*) = (J_\perp - J_*) W_L$, where $W_L$ is the sum of $L$ rung operators. Centering each rung operator by $+I/4$ gives $`\|W_L + \tfrac{L}{4} I\|_2 \le L/2`$. For $|J_\perp - J_*| \le \delta$, the scalar shift cancels in both purities and Weyl's inequality bounds the variation of each centered log-partition function by $\frac{1}{2}\beta \delta L$. Thus, given certified defects $(\varepsilon_{\mathrm{s}}, \varepsilon_\tau)$ at $(n_0, \beta_0, J_*)$, valid upper defect bounds across the entire interval $|J_\perp - J_*| \le \delta$ are

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
