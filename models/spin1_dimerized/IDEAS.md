# Candidate routes: dimerized spin-1 Heisenberg ring

Grouping adjacent spin-1 sites joined by a weak bond $1 - \delta$ into a two-spin unit cell
(dimension $3^2 = 9$) yields a translation-invariant ring of $L$ cells coupled by strong bonds
$1 + \delta$. At $L = 6$ cells ($12$ spins, Hilbert-space dimension $3^{12} = 531{,}441$), the
largest unreduced $S^z$ block has dimension $73{,}789$; exploiting two-site translation, spatial
reflection, and spin-flip symmetries reduces the largest block to dimension $11{,}596$, which can
be diagonalized densely. Within a 6-cell ($12$-spin) ring ceiling, the available seed sizes are
$n_0 = 2$ cells (requiring $2$ and $4$ cells) and $n_0 = 3$ cells (requiring $3$ and $6$ cells),
whereas $n_0 = 4$ cells would require a $16$-spin ring ($3^{16} = 43{,}046{,}721$ states).

## Candidate routes

1. **Reflected unit-cell RP identity and positive cell-transfer factorization ($X_\beta = K K^*$):**
   In the Cartesian spin-1 basis ($G_\alpha = i S^\alpha$ real and antisymmetric), every bond
   satisfies $-J\,\mathbf{S}_i\cdot\mathbf{S}_j = J \sum_\alpha G_\alpha \otimes G_\alpha$ with
   $J > 0$. Inside a two-spin cell $(a, b)$, the weak bond $o = (1 - \delta)\,\mathbf{S}_a\cdot\mathbf{S}_b$
   is real symmetric and swap-invariant ($F = P$), while the strong inter-cell bond satisfies
   $-h = (1 + \delta)\sum_\alpha C_\alpha \otimes F \overline{C_\alpha} F$ with $C_\alpha = I \otimes G_\alpha$.
   Crucially, because both incident bonds at any single site insert the exact same real matrices
   $G_\alpha$, the single-site event trace kernel $F_{\mathrm{site}}(a, u)$ is real and symmetric.
   Cutting the ring at the strong bonds and integrating out the weak intra-cell bond history $u$
   expresses the cell transfer kernel as $X_\beta(a, b) = \int F_{\mathrm{site}}(a, u) F_{\mathrm{site}}(u, b)\,d\mu_-(u) = (K K^*)(a, b)$.
   Thus $X_\beta$ is positive semidefinite (not merely self-adjoint), so all of its eigenvalues are
   nonnegative and the spatial squaring and moment-interpolation inequalities hold for **odd**
   $n_0 = 3$ cells as well as even $n_0 = 2$. **Chosen theoretical route.**
2. **Verified dense diagonalization with exact integer dihedral projectors:** In the $S^z$ product
   basis, every off-diagonal exchange matrix element of the spin-1 Heisenberg bond is equal to the
   bond coupling $J_j$ (the $\sqrt{2}$ ladder factors cancel the $1/2$), so $\mathrm{den}(\delta)\cdot H_L$
   is an exact symmetric integer matrix. For cell counts $L \in \{2, 3, 4, 6\}$, twice every cosine
   $2\cos(2\pi k r / L)$ in the real dihedral (translation + reflection) character projector is an
   exact integer. We can therefore verify projector symmetry, idempotence, and exact rank over
   $\mathbb{Z}$ on every permutation orbit, bound the distance from the floating-point orthonormal
   basis to the exact projector range, and enclose all block eigenvalues via polar-factor/Weyl
   bounds. **Chosen computational route.**
3. **Strong-dimer cluster expansion around $\delta = 1$:** A convergent low-temperature polymer
   expansion around decoupled spin-1 singlets at $\delta = 1$ could reach larger $n_0$ at smaller
   $\delta$, but deriving explicit finite-temperature remainder bounds is not needed to certify
   $\delta = 9/10$ and $\delta = 7/10$.

## Outcome

Combining Route 1 and Route 2 certified $\delta = 9/10$ at the even seed $(n_0, \beta_0) = (2, 7/2)$
($\Delta_L \ge 0.528801949327$ for even $L \ge 2$ cells) and $\delta = 7/10$ at the odd seed
$(n_0, \beta_0) = (3, 17/4)$ ($\Delta_L \ge 0.218794169839$ for even $L \ge 3$ cells), using exact
integer dihedral projectors across all 79 blocks of the 12-spin ring (`METHOD.md`). Within the
6-cell ($12$-spin) ring cap, $\delta = 3/5$ and $\delta = 1/2$ do not yet reach the admissible
region; finite-size energy-defect extrapolations place their first admissible seeds near $n_0 \approx 5$
and $n_0 \approx 7$ cells ($20\text{--}28$ spins).
