# Candidate routes: dimerized spin-1/2 Heisenberg chain

Grouping each strong bond $`(1 + \delta)\,\mathbf{S}_{2j-1}\cdot\mathbf{S}_{2j}`$ into a two-spin
unit cell produces a translation-invariant ring of $L$ cells connected by weak bonds
$`(1 - \delta)\,\mathbf{S}_{2j}\cdot\mathbf{S}_{2j+1}`$. Because a seed at $n_0$ two-spin cells
requires evaluating the thermal purity $P_\tau(2n_0, \beta_0)$ on a doubled ring of $2n_0$ cells
($4n_0$ spins), a 16-spin exact diagonalization ceiling permits even seed lengths $`n_0 \in \{2, 4\}`$
cells (while a 20-spin diagnostic ceiling reaches $n_0 = 5$ cells).

## Candidate routes

1. **Reflected unit-cell RP identity (`-h = Σ C_α ⊗ conj(C_α)` in a complex singlet/triplet basis):**
   An unreflected bond decomposition $\sum_\alpha C_\alpha \otimes C_\alpha$ ($F = I$) must commute
   with swapping two entire unit cells ($F_{\mathrm{outer}}$), which fails for the oriented
   inter-cell bond $`\mathbf{S}_{2j}\cdot\mathbf{S}_{2j+1}`$. However, combining the internal
   two-spin swap $F$ inside each cell with single-site rotation $J = 2iS^y$ on the second spin and
   transforming to the complex singlet/triplet basis $M/\sqrt{2}$ diagonalizes $F$, turning the
   antiunitary reflection $F\overline{(\cdot)}F$ into standard entrywise complex conjugation
   $-h = \sum_\alpha C_\alpha \otimes \overline{C_\alpha}$ with real diagonal on-site term $o$.
   **Chosen RP formulation.**
2. **Verified dense diagonalization in $(S^z, k, p)$ symmetry blocks:** Total magnetization $S^z$,
   one-cell translation momentum $k$, and global spin flip $p$ commute with $H_L$, decomposing the
   $2^{2n}$-dimensional Hilbert space into orthogonal symmetry blocks whose dimensions and
   multiplicities sum to $2^{2n}$. Enclosing block entries in `IntervalArithmetic` and bounding
   dense `Float64` eigenvector residuals and orthogonality defects under an a priori binary64
   rounding analysis (now the core `verified_eigen`) certifies complete spectra at $n_0 = 2$ and
   $n_0 = 4$ cells ($8$ and $16$ spins) in a few minutes.
   **Chosen route for $`\delta \in \{1/2, 3/10, 1/4\}`$.**
3. **Continuous-$`\beta`$ interval exclusion covers:** When a grid scan finds no admissible seed at
   $`n_0 \in \{2, 4\}`$ for weaker dimerizations ($`\delta \in \{1/5, 3/20, 1/10\}`$), we can upgrade
   that numerical observation to a rigorous all-$`\beta`$ exclusion theorem by combining thermal
   purity monotonicity on $(0, 1]$, a ground-state energy defect bound on $[B, \infty)$, and a
   finite cover of $[1, B]$ by interval boxes propagating lower defect bounds through the monotone
   doubling map until $z \ge 1/2$. **Chosen route for certifying the $n_0 \le 4$ grid boundary.**
4. **Operator-norm perturbation around decoupled dimers ($\delta = 1$):** At $\delta = 1$ ($w = 0$),
   the partition function factors exactly into independent dimers,
   $Z_n^{\mathrm{dimer}}(\beta) = (e^{3\beta s/4} + 3e^{-\beta s/4})^n$. Writing $H_n = H_n^{\mathrm{dimer}} + V_n$
   with $`\|V_n\| \le \frac{3}{4} n (1 - \delta)`$ and applying Weyl's inequality bounds
   $|\log Z_n - \log Z_n^{\mathrm{dimer}}| \le \frac{3}{4}\beta n (1 - \delta)$ at arbitrary $n$ in
   $O(\log n)$ operations. However, because unlinked norm bounds lose the extensive bulk free-energy
   cancellation between $\log Z_{2n}$ and $2\log Z_n$, the resulting error grows with $n$ and only
   closes seeds very close to $\delta = 1$. **Implemented as a scalable baseline and tested at 512
   spins ($\delta = 0.999999$).**
5. **Convergent linked cluster / polymer expansion around $\delta = 1$:** A Kotecký–Preiss polymer
   expansion in the weak-to-strong ratio $(1 - \delta)/(1 + \delta)$ ($1/3$ at $\delta = 1/2$,
   $7/13$ at $\delta = 3/10$) cancels the bulk free energy and bounds $\log Z_{2n} - 2\log Z_n = O(e^{-n/\xi})$
   uniformly in $n$, but requires deriving explicit low-temperature imaginary-time polymer weights
   and convergence radii.
6. **20-spin sparse Krylov diagnostics:** Using `KrylovKit` to extract the 40 lowest eigenvalues in
   large 20-spin symmetry blocks ($n_0 = 5$ cells) provides numerical insight into how the seed
   defects and energy-defect decay length $\xi_{\mathrm{eff}}$ scale with $\delta$, without
   entering the certified proof path.
7. **Coupling identity with Rayleigh–Ritz lower bounds on the doubled ring ($n_0 = 6$ cells):**
   The complete spectrum of the doubled ring is needed only for an *upper* bound on
   $Z_{2n_0}(\beta_0)$ in $\varepsilon_\tau$. With the exact coupling identity,
   ```math
   P_\tau(2n_0, \beta) = \frac{Z_{2n_0}(2\beta)}{Z_{2n_0}(\beta)^2}
     = \frac{Z_{2n_0}(2\beta)}{Z_{n_0}(\beta)^4}\, P_{\mathrm{s}}(n_0, \beta)^{-2}
     \ge \frac{Z_{2n_0}(2\beta)}{Z_{n_0}(\beta)^4}
   ```
   whenever $P_{\mathrm{s}}(n_0, \beta) \le 1$ (reflection positivity at even $n_0$), so both
   defects need only the complete verified spectrum of the ring of $n_0$ cells and *lower*
   bounds on $Z_{2n_0}(\beta_0)$ and $Z_{2n_0}(2\beta_0)$. Any set of rigorous upper bounds on low
   eigenvalues gives such lower bounds, and Rayleigh–Ritz values of approximate Lanczos vectors
   in the sparse charge and momentum blocks are rigorous upper bounds (Courant–Fischer) once
   $V^\top H V$ and $V^\top V$ carry a priori rounding bounds (new core `RitzBounds`). Omitted
   states only lower the bound, so the eigensolver is never trusted.
   *Cost of the defect:* the thermal bound is $1 - P_\tau(2n_0)\,P_{\mathrm{s}}(n_0)^2$, i.e.
   $\varepsilon_\tau^{\mathrm{coupled}} \approx \varepsilon_\tau + 2\varepsilon_{\mathrm{s}}$. At
   $n_0 = 4$ this is fatal: at $\delta = 1/4$, $\beta_0 = 29/4$ the coupled defects are
   $(0.166, 0.339)$ (computed with the core route from the complete 16-spin spectrum), which are
   not admissible. The route needs a seed length at which $\varepsilon_{\mathrm{s}}$ itself is
   small, i.e. $n_0 = 6$ cells: complete spectrum of 12 spins (4096 states, seconds) and Ritz
   bounds on the ring of 24 spins.
   *Scaling:* the 24-spin ring splits into charge (total $S^z$) sectors of at most
   $\binom{24}{12} = 2.7 \cdot 10^6$ states and $12$ momenta (plus the flip parity at $S^z = 0$),
   giving sparse blocks of dimension $0.8\text{--}2.1 \cdot 10^5$ with about $25$ entries per
   row. Bounding $40$ levels in a block takes about one minute with one thread (Lanczos,
   `KrylovKit`), so a full $\delta$ point takes under an hour. The low-temperature weight sits in
   the sectors $|S^z| \le 2$ (zero, one and two magnons, gap $\approx 0.6$ at $\delta = 1/5$).
   A first run kept $40, 40, 30, 12$ levels per block at $|S^z| = 0, 1, 2, 3$ (46 minutes per
   $\delta$); truncating it to 10 levels per block, or dropping $|S^z| = 3$, changed no defect
   in the fifth decimal near the optimal $\beta_0$. The certificates therefore keep
   $16, 16, 10$ levels at $|S^z| = 0, 1, 2$ (about 14 minutes per $\delta$, which also makes the
   `--full` replay affordable); the dropped weight enters only as a slightly larger defect.
   $n_0 = 8$ cells needs a 32-spin doubled ring ($4^{16} \approx 4.3 \cdot 10^9$ states, blocks
   of $\approx 3.7 \cdot 10^7$): out of reach of the sector enumeration in `ring_blocks`
   (capped at $2^{26}$ states) and of the memory budget. An odd $n_0 = 5$ (10 and 20 spins) would
   need $P_{\mathrm{s}}(5, \beta) \le 1$, i.e. a positive transfer operator, which this model
   lacks a proof of (`code/models/spin1_dimerized` has one for its own cell structure); it is
   only worth pursuing if $n_0 = 6$ fails. **Chosen route: it improves every certified $\delta$ and
   reaches $\delta = 1/5$.**

## Outcome

- **Exact RP:** Verified with `RingModels.verify_rp` at all six $\delta$ values in the reflected
  form with the cell swap, in the sublattice-rotated frame where the cell charge $z_a - z_b$ is
  conserved (and, as an audit, in the conjugate singlet/triplet form, `src/rp.jl`).
- **Certified positive bounds ($n_0 = 4$ cells):** Route 2, now run entirely by the shared core
  (`verified_ring_spectrum` by charge, momentum, and charge flip; `certify_gap`), certifies
  $\Delta_L \ge 0.534794$ ($\delta = 1/2, \beta_0 = 6$), $\Delta_L \ge 0.173204$
  ($\delta = 3/10, \beta_0 = 13/2$), and $\Delta_L \ge 0.058232$ ($\delta = 1/4, \beta_0 = 29/4$) for
  all even $L \ge 4$ cells (interval bound; the exact rational chain gives
  $0.534267, 0.173082, 0.058177$).
- **Rigorous all-$`\beta`$ exclusions ($n_0 = 2, 4$ cells):** Continuous-$`\beta`$ interval covers
  prove that no $\beta_0 > 0$ yields an admissible seed for this recurrence at $n_0 = 2$ (for any of
  the six $\delta$ values) or at $n_0 = 4$ for $`\delta \in \{1/5, 3/20, 1/10\}`$, establishing
  $\delta = 1/4$ as the exact certifiability threshold on rings of at most 16 spins.
- **Beyond-ED and 20-spin diagnostics:** The arbitrary-$`n`$ dimer perturbation bound certifies a
  512-spin seed at $\delta = 0.999999$ in the unit tests while quantifying why unlinked norm bounds
  fail at $\delta \le 1/2$. The 20-spin Krylov diagnostics identify a numerical odd-seed candidate
  at $(\delta, n_0, \beta_0) = (1/5, 5, 8.5)$ and confirm that the next even seeds for
  $`\delta \in \{1/5, 3/20, 1/10\}`$ require $n_0 \ge 6\text{--}8$ cells ($24\text{--}32$ spins).
- **Coupled seeds at $n_0 = 6$ cells (route 7):** the new core `RitzBounds` (sparse blocks,
  rigorous Rayleigh–Ritz upper bounds) and `coupled_seed_bounds` certify, for every even
  $L \ge 4$ cells, $\Delta_L \ge 0.657292$ ($\delta = 1/2$, $\beta_0 = 8$), $0.273433$
  ($3/10$, $65/8$), $0.168489$ ($1/4$, $17/2$), and $0.035841$ ($1/5$, $19/2$; the first
  certificate at this dimerization), up from $0.534794, 0.173204, 0.058232$. The bounds are
  $0.55, 0.34, 0.24, 0.06$ of the numerical 20-spin gaps. The coupling costs about
  $2\varepsilon_{\mathrm{s}}$ in the thermal defect: at $\delta = 1/5$ the seed is
  $(\varepsilon_{\mathrm{s}}, \varepsilon_\tau) = (0.066, 0.158)$, close to the edge of the
  admissible region. At $\delta = 3/20$ the same route fails (floating-point scan of the
  rigorous coupled defects, $\beta_0 \in [6, 15]$: $\varepsilon_\tau \ge 0.317$ throughout, with
  $\varepsilon_{\mathrm{s}} \ge 0.068$), and $\delta = 1/10$ is further out; both would need
  $n_0 = 8$ (32 spins), a sharper thermal bound than the coupling identity (an upper bound on
  $Z_{2n_0}(\beta_0)$), or an odd seed with a positivity proof.
