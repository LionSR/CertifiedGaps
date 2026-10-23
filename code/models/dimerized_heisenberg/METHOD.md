# Method: rotated cell frame, core certificates, and exclusion arguments

This document details the exact frame in which the dimerized chain is a reflection-positive core
`RingModel`, how the shared pipeline `CertifiedGaps.GapCertificates` certifies its seeds, the
continuous-$`\beta`$ exclusions, and the arbitrary-$`n`$ dimer perturbation bound. The
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
`Rational{BigInt}` with the core `is_frame_change(physical_model(δ), rotated_model(δ), W)`, and
the tests also compare the verified 8-spin spectrum with the `Float64` spectrum
(`float_ring_spectrum`) of the physical chain. All entries of $o$ and $h$ are real rationals.

**Conserved charge and flip.** The rotated exchange
$-S^xS^x + S^yS^y = -\tfrac12(S^+S^+ + S^-S^-)$ raises or lowers both spins of a bond together, and every bond joins an $a$ spin with
a $b$ spin, so the cell charge $z_a - z_b$ ($[0, 1, -1, 0]$ on
$|{\uparrow\uparrow}\rangle, |{\uparrow\downarrow}\rangle, |{\downarrow\uparrow}\rangle, |{\downarrow\downarrow}\rangle$) is
conserved. The cell spin flip $\sigma^x \otimes \sigma^x$ commutes with $o$ and $h$ and
reverses it. The lattice reflection $x \mapsto -x$ exchanges the two spins of every cell, which
also reverses the charge, so its on-site part is the cell swap composed with the spin flip,
$`P = F\,(\sigma^x \otimes \sigma^x)`$ (fixing $|{\uparrow\downarrow}\rangle$ and $|{\downarrow\uparrow}\rangle$, exchanging
$|{\uparrow\uparrow}\rangle$ and $|{\downarrow\downarrow}\rangle$): it preserves the charge, commutes with the flip, fixes $o$,
and maps $h$ to its mirror image. The core `RingModel(…; charge, flip, reflection)` checks all
of this exactly; `verified_ring_spectrum` then splits every ring by charge and cell momentum,
diagonalizes the sectors $q$ and $-q$ once and the momenta $\pm k$ once, and, the model being
real, builds every block at $0 < k < \pi$ as a real symmetric matrix (smaller rounding radii).

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

## Coupled seeds at $n_0 = 6$ cells (core `RitzBounds`)

The doubled ring of $2n_0 = 12$ cells has $2^{24} \approx 1.7 \cdot 10^7$ states, too many for a
complete verified spectrum. The criterion needs
$\varepsilon_{\mathrm{s}} \ge 1 - P_{\mathrm{s}}(6, \beta_0)$
and $\varepsilon_\tau \ge 1 - P_\tau(12, \beta_0)$, and by the coupling identity

```math
P_\tau(2n_0, \beta) = \frac{Z_{2n_0}(2\beta)}{Z_{2n_0}(\beta)^2}
  = \frac{Z_{2n_0}(2\beta)}{Z_{n_0}(\beta)^4}\, P_{\mathrm{s}}(n_0, \beta)^{-2}
  \ge \frac{Z_{2n_0}(2\beta)}{Z_{n_0}(\beta)^4},
```

since $P_{\mathrm{s}}(n_0, \beta) = \mathrm{Tr}\,X^{2n_0}/(\mathrm{Tr}\,X^{n_0})^2 \le 1$ for the
self-adjoint transfer operator $X$ of the reflection-positive ring and even $n_0$. So the core
`coupled_seed_bounds` takes

```math
\varepsilon_{\mathrm{s}} = 1 - \frac{\underline{Z_{12}}(\beta_0)}{\overline{Z_6}(\beta_0)^2}, \qquad
\varepsilon_\tau = 1 - \frac{\underline{Z_{12}}(2\beta_0)}{\overline{Z_6}(\beta_0)^4},
```

with $\overline{Z_6}$ from the complete verified spectrum of the ring of 6 cells (4096 states,
`verified_ring_spectrum`) and lower bounds $\underline{Z_{12}}$ from rigorous upper bounds on the
low levels of the 24-spin ring:

* `ring_blocks(model, 12; sparse = true, charges = 0:2)` builds the charge (total $S^z$),
  momentum and flip-parity blocks of the rotated model exactly as for the dense certificates
  (interval-enclosed entries, real basis at $0 < k < \pi$, exact multiplicities, an exact
  Frobenius bound on the distance of the `Float64` center to the true block), but stored sparse
  ($0.8\text{--}2.1 \cdot 10^5$ rows, about 25 entries per row). The sectors $\pm q$ share a
  block through the flip; $|S^z| \ge 3$ is omitted.
* In each block KrylovKit's Lanczos (deterministic start vector) gives approximate eigenvectors
  $V$ of the $k = 16, 16, 10$ lowest levels at $|S^z| = 0, 1, 2$. `ritz_upper_bounds`
  turns them into rigorous upper bounds: by Courant–Fischer the Ritz values of $\mathrm{span}\,V$
  bound the $k$ lowest eigenvalues from above; $V^\top M V$ and $V^\top V$ are computed in `Float64`
  with a priori $\gamma_n$ bounds (sparse product with $\gamma_{p+2}\,\| |M| \|_\infty$, GEMMs
  with $\gamma_{N+2}$); the $k \times k$ matrix is enclosed by `verified_eigen`; and the pencil
  $(V^\top M V, V^\top V)$ is reduced to it by a Gershgorin shift and Ostrowski's theorem (the
  derivation is the docstring of `ritz_upper_bounds`). The block's construction error is added.
* `log_partition_lower` sums $m_B e^{-\beta u_i}$ in interval arithmetic. Every omitted state
  only lowers the bound; the Lanczos vectors are never trusted.

The defects are then raised by $10^{-9}$ and rounded up to that grid (`padded_seed`), so a
replay whose Lanczos vectors differ in the last bits (another BLAS) still checks, and
`certify_point` writes a certificate with `seed.route = "coupled"`, which also pins the SHA-256 of
the stored bounds (`results/ritz_bounds_12_cells_delta_*.json.gz`). The seed theorem covers every
even $L \ge 6$ cells; the remaining even length $L = 4$ is covered by the verified gap of the
4-cell ring, which exceeds the bound. The price of the coupling is
$\varepsilon_\tau \approx \varepsilon_\tau^{\mathrm{true}} + 2\varepsilon_{\mathrm{s}}$,
which rules out
$n_0 = 4$ (at $\delta = 1/4$: $(0.166, 0.339)$, not admissible) but is small at $n_0 = 6$.

## Rigorous all-$`\beta`$ exclusions (core `Exclusions`)

To prove that no $\beta_0 > 0$ satisfies the doubling recurrence at a given $(n_0, \delta)$,
`certify.jl` calls the core `Exclusions.exclude_all_beta` on the verified spectra of the rings
$n_0$ and $2n_0$:
1. **Tail $[B, \infty)$:** with
   $`d_E = \underline{E_0(2n_0)} - 2\,\overline{E_0(n_0)} > 0`$ from `ground_energy_enclosure`
   and $`F_L(\beta) = \log \mathrm{Tr}\, e^{-\beta (H_L - a_L)}`$ ($a_L$ the lower end of
   $E_0(L)$, so $F_L$ decreases), $\log P_{\mathrm{s}}(n_0, \beta) \le F_{2n_0}(B) - B d_E$ for all
   $\beta \ge B$, with $B$ the smallest power of two whose lower spatial defect is rejected.
2. **Windows covering $[0, B]$:** $[0, B]$ is bisected into dyadic windows. On each, monotonicity of
   $F_L$ gives a rational lower bound on $1 - P_{\mathrm{s}}$, and
   $\frac{d}{d\beta}\log P_\tau(2n_0, \beta) = 2(\langle H \rangle_\beta - \langle H \rangle_{2\beta}) \ge 0$ gives
   $1 - P_\tau(2n_0, \beta) \ge 1 - P_\tau(2n_0, \beta_{\mathrm{hi}})$; these lower defects are pushed
   through the monotone doubling map (`Exclusions.rejection`, rounding down) until an argument
   $z$ reaches $1/2$, after which the recurrence cannot converge.

The summary stores every window with its lower bounds and rejection witness and the tail; the fast
check (`check_exclusion`) reruns the recurrence on them and checks that the windows tile
$[0, B]$, and `check.jl --full` recomputes every bound from fresh spectra.

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
