# Route memo: spin-1 chain with single-ion anisotropy $D$

We study the periodic spin-1 chain $`H_L = \sum_{j=1}^L \mathbf{S}_j \cdot \mathbf{S}_{j+1} + D \sum_{j=1}^L (S_j^z)^2`$
with single-site unit cells in the Cartesian reflection-positive basis. This memo records the
candidate routes evaluated prior to Stage 1 (the verified exact-diagonalization baseline) and
Stage 2 (the transfer-squaring and Weyl-transport extensions), followed by their outcomes.

## Candidate routes

1. **Closed forms:** Because the nearest-neighbor Heisenberg exchange terms commute neither with the
   single-ion anisotropy $D(S_j^z)^2$ nor with adjacent bonds, no closed-form or free-fermion
   partition function is available at finite $D$.
2. **Cluster / low-temperature polymer expansion around $D = \infty$:** The product state
   $|0\rangle^{\otimes L}$ ($S_j^z = 0$) is the exact ground state at $D = \infty$, with local
   excitation energy $D$, but bounding $\log Z_{2n} - 2\log Z_n$ uniformly in $n$ near $D \approx 3.2$
   against virtual $m = \pm 1$ pair creation requires a convergent low-temperature polymer
   expansion with explicit quantitative remainders.
3. **Symmetry-resolved verified dense exact diagonalization (Stage-1 baseline):** Conservation of
   total $S^z$, translation momentum, spin inversion at $S^z = 0$, and conjugate-momentum pairing
   reduce the $3^{12} = 531{,}441$-dimensional Hilbert space at $L = 12$ into 92 blocks of dimension
   at most $5{,}798$. Every block can be diagonalized densely and certified via interval matrix
   construction, IEEE binary64 residual and orthogonality norm bounds, and polar-factor/Weyl
   eigenvalue enclosures, yielding rigorous seeds at $n_0 = 6$ ($2n_0 = 12$). Going directly to
   $n_0 = 8$ ($2n_0 = 16$) or even $L = 14$ (largest block $41{,}835$) is prohibitive for dense
   diagonalization. **Chosen Stage-1 route.**
4. **Spatial transfer squaring and purity coupling without the $4n_0$ ring (Stage-2 extension):**
   For even $n$, $\sigma_n = X_\beta^n / Z_n(\beta)$ is a positive trace-one density operator, and
   squaring bounds the spatial defect at $2n$ by $\varepsilon_{\mathrm{s}}(2n, \beta) \le f(\varepsilon_{\mathrm{s}}(n, \beta)) = \varepsilon_{\mathrm{s}}^2 / [2(1 - \varepsilon_{\mathrm{s}})^2]$.
   Meanwhile, the exact partition-function identity
   ```math
   P_\tau(4n, \beta) = \frac{P_{\mathrm{s}}(2n, 2\beta)\,P_\tau(2n, \beta)^2}{P_{\mathrm{s}}(2n, \beta)^2} \ge P_{\mathrm{s}}(2n, 2\beta)\,P_\tau(2n, \beta)^2
   ```
   bounds the thermal purity at $4n = 24$ using only $\log Z_6$ and $\log Z_{12}$ evaluated at
   $\beta$ and $2\beta$ (equivalent to Theorem 2.24(v) of the survey). This lifts an $n_0 = 6$
   spectral calculation to an $n_0 = 12$ (or $24, 48, 96$) seed in $O(\log(n_0/6))$ interval
   operations without ever diagonalizing $L = 24$. **Chosen Stage-2 route.**
5. **Centered Weyl perturbation transport in $D$ (Stage-2 extension):** Because
   $V = \sum_{j=1}^L (S_j^z)^2$ satisfies $0 \le V \le L I$, centering $V$ by $\frac{L}{2}I$ gives
   $`\|\delta D (V - \frac{L}{2}I)\|_2 \le \frac{1}{2}L|\delta D|`$. By Weyl's inequality, shifting
   $`D^* \mapsto D^* + \delta D`$ shifts every eigenvalue of $`H_L(D) - \frac{L\delta D}{2}I`$ by at
   most $`\frac{1}{2}L|\delta D|`$, and the constant shift cancels in both purities. This transports
   verified $`\log Z`$ enclosures across narrow intervals of $`D`$ without re-diagonalizing $`H_{12}`$.
   **Chosen Stage-2 route.**

## Outcome

- **Stage-1 ED seeds ($n_0 = 6$):** Complete verified $L = 6, 12$ spectra certify $D = 4$
  ($\Delta_L \ge 0.694746580$ at $\beta_0 = 9/4$) and $D = 1603/500 = 3.206$
  ($\Delta_L \ge 0.002752708$ at $\beta_0 = 381/125$) for all even $L \ge 6$. On the exploratory
  grid, $D \le 3.205$ has no admissible even seed at $n_0 \le 6$.
- **Stage-2 transfer squaring ($n_0 = 12, \dots, 96$):** At $D = 4$, lifting from $L = 6, 12$ to a
  seed at $(n_0, \beta_0) = (12, 29/10)$ without constructing $L = 24$ raises the certified uniform
  gap bound from $0.694746580$ to **$0.899422073$** for all even $L \ge 12$, and verified spectra at
  $L = 6, 8, 10$ bridge the remaining even lengths down to $L \ge 6$. Iterating further certifies
  $n_0 = 24, 48, 96$ as well. At $D = 3.206$, where the $n_0 = 6$ thermal defect is already near
  threshold, spatial lifting without a new thermal input loses too much thermal purity to close.
- **Stage-2 Weyl transport:** Centered Weyl perturbation around $D^* = 3.206$ certifies the
  continuous parameter interval $D \in [3.2059998, 3.206]$ ($\Delta_L \ge 0.001821437$ for all even
  $L \ge 6$), and similarly certifies $D \in [3.9979773, 4]$ around $D^* = 4$.
- **Consolidation onto the shared core:** Stage 1 now runs the core pipeline
  (`RingModel` in the $S^z$ frame with the reflected RP form and spin flip, `certify_gap`,
  `verified_eigen`, `write_spectrum`), and stage 2 reads the stored core spectra. The core
  enclosure is tighter than the former local one, so every defect decreased and every bound
  increased; the numbers above are the regenerated ones (see `README.md`, `## History`).
