# Method: rotated cell frame, core certificates, and exclusion arguments

This document details the exact frame in which the dimerized chain is a reflection-positive core
`RingModel`, how the shared pipeline `CertifiedGaps.GapCertificates` certifies its seeds, the
continuous-$`\beta`$ exclusion covers, and the arbitrary-$`n`$ dimer perturbation bound. The
model-specific code is `src/DimerizedHeisenberg.jl`; everything else is the core.

## The rotated cell frame

Group each strong bond into a two-spin unit cell $(a, b)$, and set $s = 1 + \delta$,
$w = 1 - \delta$. The physical cell term and inter-cell bond (second spin of a cell, first spin
of the next, on $\mathrm{cell} \otimes \mathrm{cell}$) are
```math
o_{\mathrm{phys}} = s\, \mathbf{S}_a \cdot \mathbf{S}_b, \qquad
h_{\mathrm{phys}} = w\, (I \otimes \mathbf{S}) \cdot (\mathbf{S} \otimes I).
```
Rotate the second spin of every cell by
$`J = 2iS^y = \begin{psmallmatrix} 0 & 1 \\ -1 & 0 \end{psmallmatrix}`$, a real orthogonal matrix with $J^\top S^x J = -S^x$, $J^\top S^z J = -S^z$,
and $J^\top (iS^y) J = iS^y$. With $W = I \otimes J$ and $U = W^{\otimes L}$, every bond of the
ring joins a rotated and an unrotated spin, so
```math
U^\top H_L U = \sum_j h_{j,j+1} + \sum_j o_j, \qquad
o = s\,(-S^xS^x + S^yS^y - S^zS^z), \qquad
h = w\, I \otimes (-S^xS^x + S^yS^y - S^zS^z) \otimes I.
```
Because $U$ is a product of cell unitaries, the two local identities
$W^\top o_{\mathrm{phys}} W = o$ and $(W \otimes W)^\top h_{\mathrm{phys}} (W \otimes W) = h$
imply unitary equivalence for every $L$; `frame_identity` checks both exactly in
`Rational{BigInt}`, and the tests also compare the verified 8-spin spectrum with an independent
`Float64` solver of the physical chain. All entries of $o$ and $h$ are real rationals.

**Conserved charge and flip.** The rotated exchange
$-S^xS^x + S^yS^y = -\tfrac12(S^+S^+ + S^-S^-)$ raises or lowers both spins of a bond together, and every bond joins an $a$ spin with
a $b$ spin, so the cell charge $z_a - z_b$ ($[0, 1, -1, 0]$ on
$|{\uparrow\uparrow}\rangle, |{\uparrow\downarrow}\rangle, |{\downarrow\uparrow}\rangle, |{\downarrow\downarrow}\rangle$) is
conserved. The cell spin flip $\sigma^x \otimes \sigma^x$ commutes with $o$ and $h$ and
reverses it. The core `RingModel(…; charge, flip)` checks both exactly; `verified_ring_spectrum`
then splits every ring by charge and cell momentum and diagonalizes the sectors $q$ and $-q$
once.

## Reflected reflection positivity

Let $F$ be the swap of the two spins of a cell (real, symmetric, $F^2 = I$, $FoF = o$). With
$`P_\alpha \in \{S^x, iS^y, S^z\}`$ (real),
```math
-h = w \sum_{\alpha} C_\alpha \otimes F\, \overline{C_\alpha}\, F, \qquad
C_\alpha = I \otimes P_\alpha, \qquad F (I \otimes P_\alpha) F = P_\alpha \otimes I,
```
which is the reflected form of `notes/thermal-mpo-gap-scope.tex` (Appendix E) and of the Lean
structure `RPRing` (`lean/CertifiedGaps/Models/DimerizedHeisenberg.lean` uses exactly this
frame). `RingModels.verify_rp` checks it exactly, and every certificate stores the factors, the
weights $w$, and the reflection $F$. Without $F$ there is no such factorization: the core's exact
LDL$`^\dagger`$ factorization of the realigned bond fails in this basis, and
`src/reflection_conventions.jl` shows that no common single-cell basis change helps, since
$`\|h - F_{\mathrm{outer}} h F_{\mathrm{outer}}\|_F^2 = \tfrac38 (1 - \delta)^2 > 0`$ for the swap
$`F_{\mathrm{outer}}`$ of two whole cells is basis independent.

`src/rp.jl` (`results/rp_audit.json`) records an equivalent second frame: the singlet/triplet
basis $M/\sqrt{2}$ with
```math
M = \begin{pmatrix} 1 & 1 & 0 & 0 \\ 0 & 0 & 1 & i \\ 0 & 0 & 1 & -i \\ 1 & -1 & 0 & 0
\end{pmatrix}, \qquad M^\dagger M = 2I, \qquad \tfrac12 M M^\top = F,
```
diagonalizes $F$ and turns the reflected form into the conjugate form
$`-h = \sum_\alpha C'_\alpha \otimes \overline{C'_\alpha}`$ with complex $C'_\alpha$ and
$o = \tfrac{s}{4}\mathrm{diag}(-3, 1, 1, 1)$. Both forms are verified exactly; the
certificates use the real rotated frame because it carries the conserved charge.

## Core certificates

For each seed $(n_0, \beta_0)$, `certify.jl` computes the verified spectra of the rings of $n_0$
and $2n_0$ cells once with `VerifiedSpectra.verified_ring_spectrum` and passes them to
`GapCertificates.seed_bounds` and `certify_gap`:

* every charge and momentum block is built exactly from the rational model (phases and orbit
  ratios enclosed in interval arithmetic); its binary64 center and an exact Frobenius bound on
  the construction error go to `verified_eigen`, which bounds the residual and orthogonality
  defects of the LAPACK eigenvectors under the IEEE binary64 rounding model
  ($\gamma_{N+8}$ for real blocks, $\gamma_{2N+8}$ for complex blocks evaluated with four real
  GEMMs, gradual-underflow allowances, Frobenius norms) and returns the smaller of the
  Ostrowski and polar-factor Weyl radii;
* $\log Z_{n_0}(\beta_0)$, $\log Z_{2n_0}(\beta_0)$, and $\log Z_{2n_0}(2\beta_0)$ are enclosed by
  `log_partition_enclosure`, and the defects
  $1 - P_{\mathrm{s}}(n_0, \beta_0) = 1 - e^{\log Z_{2n_0}(\beta_0) - 2\log Z_{n_0}(\beta_0)}$ and
  $1 - P_\tau(2n_0, \beta_0) = 1 - e^{\log Z_{2n_0}(2\beta_0) - 2\log Z_{2n_0}(\beta_0)}$ are rounded
  up to the grid $10^{-15}$ to give $(\varepsilon_{\mathrm{s}}, \varepsilon_\tau)$;
* `certified_gap_bound` screens the seed, and `rational_gap_chain` gives the exact rational
  chain whose final bound is the certificate's strict `gap_lower_bound`.

`GapCertificates.check_certificate` re-derives the RP identity and the chain from the stored
seed (and, with `spectra = true`, the defects from fresh spectra). The 16-spin gap interval is
`spectral_gap_enclosure` of the 8-cell spectrum, rounded outward to $10^{-6}$.

## Rigorous all-$`\beta`$ exclusion covers (`negative_cover`)

To prove that no $\beta_0 > 0$ satisfies the doubling recurrence at a given $(n_0, \delta)$,
`negative_cover` partitions $(0, \infty)$ into three regimes, using only core enclosures of the
verified spectra of the rings $n_0$ and $2n_0$:
1. **High temperature $(0, 1]$:**
   $\frac{d}{d\beta}\log P_\tau(2n_0, \beta) = 2(\langle H \rangle_\beta - \langle H \rangle_{2\beta}) \ge 0$, so $P_\tau(2n_0, \beta)$ is
   nondecreasing; $P_\tau(2n_0, 1) \le 1/2$ (`log_partition_enclosure` at $\beta = 1$) gives
   $\varepsilon_\tau \ge 1/2$, hence $z \ge 1/2$, for all $\beta_0 \le 1$.
2. **Low temperature $[B, \infty)$:** with
   $`d_E = \underline{E_0(2n_0)} - 2\,\overline{E_0(n_0)} > 0`$ from `ground_energy_enclosure`,
   $P_{\mathrm{s}}(n_0, \beta) \le 2^{4n_0} e^{-\beta d_E} \le 1/2$ for all $\beta \ge B$, with $B$
   the smallest power of two $\ge 32$ that works.
3. **Intermediate $[1, B]$:** $[1, B]$ is bisected into closed boxes. On each box, interval
   evaluation of `shifted_log_partition` (one exact shift per ring, so the box width does not
   multiply the ground energies) gives rational lower bounds on $\varepsilon_{\mathrm{s}}(\beta)$
   and $\varepsilon_\tau(\beta)$; these are propagated through the monotone map
   $(\varepsilon_{\mathrm{s}}, \varepsilon_\tau) \mapsto (\varepsilon_{\mathrm{s}}', \varepsilon_\tau')$,
   rounding down to $10^{-15}$ at each step, until $z \ge 1/2$.

The summary stores every box with its lower bounds; the fast check reruns the recurrence on them
and checks that the boxes tile $[1, B]$, and `check.jl --full` recomputes every bound from fresh
spectra.

## Arbitrary-$`n`$ dimer perturbation bound (`dimer_bounds`)

At $\delta = 1$ ($w = 0$), the $n$-cell ring consists of $`n`$ decoupled strong dimers with exact
partition function $Z_n^{\mathrm{dimer}}(\beta) = (e^{3\beta s/4} + 3e^{-\beta s/4})^n$. For
$w = 1 - \delta > 0$, writing $H_n = H_n^{\mathrm{dimer}} + V_n$ gives
$`\|V_n\|_2 \le v_n := \frac{3}{4} n w`$, so Weyl's inequality implies
$|\log Z_n(\beta) - \log Z_n^{\mathrm{dimer}}(\beta)| \le \beta v_n$ and hence
```math
P_{\mathrm{s}}(n, \beta) \ge e^{-4\beta v_n}, \qquad
P_\tau(2n, \beta) \ge \left(\frac{1 + 3x^2}{(1 + 3x)^2}\right)^{2n} e^{-8\beta v_n}, \qquad x = e^{-\beta s}.
```
Evaluating these bounds requires only $O(\log n)$ interval operations independently of the
$2^{2n}$ Hilbert-space dimension (certifying $n_0 = 128$ cells / $512$ spins at
$\delta = 0.999999, \beta_0 = 8$ in `test/runtests.jl`), but because the extensive error
$4\beta v_n = 3\beta n (1 - \delta)$ does not cancel between $Z_{2n}$ and $Z_n^2$, reaching
moderate $\delta$ at large $`n`$ would require a linked cluster expansion.
