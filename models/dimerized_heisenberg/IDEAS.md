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
   dense `Float64` eigenvector residuals and orthogonality defects via a two-pass Cauchy–Schwarz
   $\gamma_{4d+16}$ rounding analysis certifies complete spectra at $n_0 = 2$ and $n_0 = 4$ cells
   ($8$ and $16$ spins) in a few minutes. **Chosen route for $`\delta \in \{1/2, 3/10, 1/4\}`$.**
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

## Outcome

- **Exact RP:** Verified in `Complex{Rational{BigInt}}` across all six $\delta$ values in the
  reflected singlet/triplet cell basis (`src/rp.jl`, `RingModels.verify_rp`).
- **Certified positive bounds ($n_0 = 4$ cells):** Complete verified 8- and 16-spin spectra
  combined with the two-pass Cauchy–Schwarz rounding bound certify $\Delta_L \ge 0.534756$
  ($\delta = 1/2, \beta_0 = 6$), $\Delta_L \ge 0.173197$ ($\delta = 3/10, \beta_0 = 13/2$), and
  $\Delta_L \ge 0.058225$ ($\delta = 1/4, \beta_0 = 29/4$) for all even $L \ge 4$ cells.
- **Rigorous all-$`\beta`$ exclusions ($n_0 = 2, 4$ cells):** Continuous-$`\beta`$ interval covers
  prove that no $\beta_0 > 0$ yields an admissible seed for this recurrence at $n_0 = 2$ (for any of
  the six $\delta$ values) or at $n_0 = 4$ for $`\delta \in \{1/5, 3/20, 1/10\}`$, establishing
  $\delta = 1/4$ as the exact certifiability threshold on rings of at most 16 spins.
- **Beyond-ED and 20-spin diagnostics:** The arbitrary-$`n`$ dimer perturbation bound certifies a
  512-spin seed at $\delta = 0.999999$ in the unit tests while quantifying why unlinked norm bounds
  fail at $\delta \le 1/2$. The 20-spin Krylov diagnostics identify a numerical odd-seed candidate
  at $(\delta, n_0, \beta_0) = (1/5, 5, 8.5)$ and confirm that the next even seeds for
  $`\delta \in \{1/5, 3/20, 1/10\}`$ require $n_0 \ge 6\text{--}8$ cells ($24\text{--}32$ spins).
