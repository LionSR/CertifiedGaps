# Spin-2 Heisenberg Chain: Route Analysis and Thermal Starting-Point Audit

## Goal and candidate routes

The spin-2 antiferromagnetic Heisenberg chain $H_L = \sum_{j=1}^L \mathbf{S}_j \cdot \mathbf{S}_{j+1}$
($\mathbf{S}_j^2 = 6I$) has a small bulk Haldane gap ($\Delta \approx 0.08917$) and a correlation length
$\xi \approx 49.5$ sites. Consequently, an admissible seed $(n_0, \beta_0)$ for the spacetime bootstrap
lies at large ring lengths and inverse temperatures (for instance, the conditional one-trace target in
`METHOD.md` uses $n_0 = 256$ and $\beta_0 = 64$, requiring $\widehat{Z}_{256}(64) \le 1.08$ relative to a
certified per-site energy shift $\kappa = 4.7612463$).

While the variational lower bound on the partition function at $L = 1024$ is already certified via an
$\mathrm{SU}(2)$-invariant MPS trial state ($D = 295$, proving $E_0(L)/L < -4.7612463$), certifying an
*upper* bound on $Z_L(\beta)$ at $L = 256, \beta = 64$ requires controlling a thermal MPO across many
imaginary-time steps. The current 256-site quadratic Taylor-majorant route starts near
$\beta \approx 10^{-7}$ and requires $28$ squarings. Before running further large-ring contractions, we
compare the candidate strategies and separate the error budget of the initial thermal approximation from
the error budget of MPO bond compression.

| Route | Mechanism and advantage | Scaling obstruction and decision |
|---|---|---|
| **1. Local Trotter gates + final-trace certification + charge blocks** | Uses the Araki–Lieb–Thirring trace inequality $Z_L(\beta) \le \operatorname{Tr}[(e^{-\beta A/m} e^{-\beta B/m})^m]$ to start at moderate step sizes $\beta/m$ without an extensive Taylor remainder, combined with signed adjoint-trace error control and $U(1)$ charge blocks. | Product-formula splitting overestimate and compressed-rank growth compete as $L$ and $\beta$ increase. **Selected as primary route**, subject to benchmarking splitting vs. compression overhead on small exact references first. |
| **2. Higher-degree Taylor / sum-of-squares polynomial majorant** | Replaces the degree-2 polynomial $Q_2(K)$ with a degree-$d$ polynomial (factored into a rational sum of squares of degree-$d/2$ polynomials) to eliminate early near-identity squarings. | While the scalar remainder drops rapidly with $d$, the exact MPO bond dimension of $K^{d/2}$ grows exponentially with $d$, shifting the bottleneck from squaring count to initial MPO construction and compression. Retained as a comparison baseline (`d = 4`). |
| **3. Correlated purity-ratio bounds or transfer-operator contraction** | Directly bounding $\log P_{\mathrm{s}}(n, \beta) = \log Z_{2n}(\beta) - 2\log Z_n(\beta)$ or contracting the spatial transfer operator $X_\beta$ would exploit cancellation of the extensive free energy $L f(\beta)$. | Independent interval bounds on $\log Z_{2n}$ and $\log Z_n$ add rather than cancel errors, and $P_\tau(2n, \beta)$ compares two different inverse temperatures ($\beta$ and $2\beta$) whose free-energy densities differ. Requires a new analytic theorem before implementation. |
| **4. Continuing the degree-2 polynomial flow with symmetry reduction alone** | Reuses the existing verified squaring pipeline with $U(1)$ charge conservation. | Reduces per-step runtime and witness size, but leaves $28$ squaring steps across which Hilbert–Schmidt or operator-norm residual bounds must remain below $O(1)$. |

## Audit of the thermal starting point

To quantify the trade-off between Route 1 (local Trotter product upper bounds) and Route 2 (polynomial
Taylor majorants), `scripts/polynomial/audit_thermal_start.jl` evaluates two exact/Arb benchmarks.

### 1. Why the degree-2 polynomial route starts at $\beta \approx 2.4 \times 10^{-7}$

The polynomial majorant route shifts $K = H_L + L I$ so that $-5L I \preceq K \preceq 5L I$. With
$m = 2^{\text{power}}$ squarings, step parameter $\delta = \beta/(2m)$ (so the initial trace temperature
is $2\delta = \beta/m$), and polynomial degree $d$, `CertifiedGaps.CertifiedMPO.PolynomialBounds` uses

$$x = 5 L \delta, \qquad r \ge \frac{e^x x^{d+1}}{(d+1)!}, \qquad Q_d(K) = \sum_{j=0}^d \frac{(-\delta K)^j}{j!} + r I \succeq e^{-\delta K} \succ 0,$$

for which the exact polynomial trace inflates the true Gibbs trace $Z_L(\beta) = e^{\beta L} \operatorname{Tr}(e^{-2\delta K})^m$
by at most $\log(\text{inflation}) \le 4 m r e^x$.

For $d = 2$, $4m r e^x \propto (\beta L)^3 / m^2$ when $x$ is small, so the extensive spectral norm
$\|K\| \le 5L$ forces $m = 2^{28}$ at $L = 256, \beta = 64$. Using 384-bit Arb arithmetic via
`taylor_seed(256, 64, p, d)`, the minimum number of squarings $p$ needed to keep the scalar trace
inflation factor below $1001/1000$ at $L = 256, \beta = 64$ is:

| Polynomial degree $d$ | Subsequent squarings $p$ | Initial trace $\beta/2^p$ |
|---:|---:|---:|
| $2$ | $28$ | $1 / 4{,}194{,}304$ |
| $4$ | $21$ | $1 / 32{,}768$ |
| $8$ | $17$ | $1 / 2{,}048$ |
| $12$ | $16$ | $1 / 1{,}024$ |
| $16$ | $15$ | $1 / 512$ |
| $32$ | $13$ | $1 / 128$ |
| $64$ | $12$ | $1 / 64$ |
| $128$ | $11$ | $1 / 32$ |
| $256$ | $10$ | $1 / 16$ |
| $512$ | $9$ | $1 / 8$ |

Raising $d$ removes up to $19$ squarings at the scalar level, but constructing and certifying a compressed
MPO for $Q_d(K)$ at high degree is itself prohibitive because the exact bond dimension of $K^k$ grows
exponentially with $k$.

### 2. Splitting vs. compression overhead in the local-gate Araki–Lieb–Thirring bound

Splitting $H_L = A + B$ into commuting even-bond and odd-bond layers, the Araki–Lieb–Thirring inequality
gives the rigorous upper bound

$$Z_L(\beta) \le \operatorname{Tr}\!\left[(e^{-\beta A/m} e^{-\beta B/m})^m\right]$$

for every positive integer $m$, with no separate Trotter remainder term. However, the product-formula
trace itself overestimates $Z_L(\beta)$ by an amount governed by $\beta/m$.

On the 4-site ring ($625$ states), the exact eigenvalues of $H_4$ are known analytically from angular-momentum
addition:

$$E(u, v, s) = \frac{s(s+1) - u(u+1) - v(v+1)}{2}, \qquad \text{multiplicity } 2s+1,$$

for $u, v \in \{0, 1, 2, 3, 4\}$ and $s \in \{|u-v|, \dots, u+v\}$. At $\beta = 4$ and per-site energy
shift $\kappa = 238061/50000 = 4.76122$, the exact shifted Gibbs trace $\widehat{Z}_4(4)$ is enclosed in
Arb arithmetic around $48.13451956631644$. Comparing `Spin2Gap.ThermalBounds.trotter_upper(4, 4, m; shift=238061//50000)` against the exact Gibbs trace gives:

| Trotter steps $m$ | Gate step $\beta/m$ | Trotter upper bound / exact Gibbs trace |
|---:|---:|---:|
| $64$ | $1/16$ | $1.054501559780$ |
| $128$ | $1/32$ | $1.013417771778$ |
| $256$ | $1/64$ | $1.003341576520$ |

In `METHOD.md`, the 5-doubling adjoint-trace experiment at $L = 4, \beta = 4$ ($m = 64$, starting at
$\beta/m = 1/16$) certified the compressed Trotter trace to within $0.03011\%$ of the *uncompressed
Trotter trace* (`results/adjoint-trace-compressed-n4.json`), but that $0.03011\%$ compression overhead
sits on top of a $5.45\%$ Trotter splitting overestimate relative to the true Gibbs trace. Both error
budgets—the product-formula overestimate as a function of $\beta/m$ and the certified MPO compression
overhead—must be small simultaneously.

## Reproduction

Both audit tables and all acceptance inequalities can be replayed from the repository root in under a
second (excluding Julia startup):

```sh
OPENBLAS_NUM_THREADS=1 julia --startup-file=no --project=models/spin2 \
  models/spin2/scripts/polynomial/audit_thermal_start.jl /tmp/spin2-thermal-start-audit.json
```

Omitting the output path prints the audit summary to stdout without writing a file. All checks use exact
integer/rational arithmetic and 384-bit Nemo/Arb ball enclosures.

## Outcome

The thermal starting-point audit (`scripts/polynomial/audit_thermal_start.jl`, 0.43 s runtime on
Julia 1.12) established two conclusions that guide further spin-2 work:
1. Increasing the initial Taylor polynomial degree $d$ reduces the number of subsequent squarings
   only by shifting the cost into an exponentially growing initial MPO construction; the quartic
   ($d = 4$) sum-of-squares seed is the practical limit of that route.
2. The local-gate Araki–Lieb–Thirring route avoids the $O(\beta L)$ Taylor step restriction, but its
   Trotter splitting overestimate ($5.45\%$ at $\beta/m = 1/16$ on $L = 4$) must be tracked
   alongside the MPO compression error.

A certified finite-size variational bound ($E_0(L)/L < -4.7612463$ at $L \in \{128, 256, 512, 1024\}$)
and the certified MPO compression toolkit (`CertifiedGaps.CertifiedMPO`) are in place, while a
low-temperature thermal trace upper bound at $L = 256, \beta = 64$—and therefore a uniform spin-2 gap
certificate—remains open.
