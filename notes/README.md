# LaTeX Notes

This directory contains the REVTeX 4.2 notes that develop the spacetime bootstrap, its scope across
reflection-positive spin chains, and the diagrammatic certification of the spin-1 Haldane gap seed. Each
`.tex` source file has its compiled `.pdf` alongside it. The main document `thermal-mpo-gap-note.tex`
also inputs `thermal-mpo-gap-scope.tex` (Appendices D–H) and `thermal-mpo-gap-certificate.tex`
(Appendices I–M), which can be compiled either as part of the main note or as standalone documents.

| note | PDF | content |
|---|---|---|
| [`thermal-mpo-gap-note.tex`](thermal-mpo-gap-note.tex) | [PDF](thermal-mpo-gap-note.pdf) | *A Diagrammatic Proof of the Haldane Gap via Spacetime Tensor Networks*. Main text: proves a uniform spectral gap $\Delta_L > (\log 20)/784 \approx 0.003821$ for the periodic spin-1 antiferromagnetic Heisenberg chain for all even $L \ge 2304$, using an exact continuous-time thermal MPO and a scale-doubling iteration on spatial and thermal purities. Appendices A–C treat open chains with equal endpoint fields $h = 3/5$, proving $\Delta_n > (\log 10)/392$ for odd $n \ge 1921$. |
| [`thermal-mpo-gap-scope.tex`](thermal-mpo-gap-scope.tex) | [PDF](thermal-mpo-gap-scope.pdf) | *Scope of the Spacetime Bootstrap* (Appendices D–H). Characterizes which Hamiltonians the method covers: nearest-neighbor bond reflection positivity (including multi-site unit cells with an internal reflection $F$) makes the spatial transfer operator $X_\beta$ self-adjoint, so bounds on the two purities at one finite scale imply a uniform spectral gap. |
| [`thermal-mpo-gap-certificate.tex`](thermal-mpo-gap-certificate.tex) | [PDF](thermal-mpo-gap-certificate.pdf) | *How the Initial Estimate Is Certified* (Appendices I–M). Presents the finite seed certificate at $(n_0, \beta_0) = (2304, 784)$ in tensor-network diagrams: twisted partition functions, symmetry-sector projectors, Chebyshev polynomial filters, a matrix-product trial vector at $n = 36$, and six unequal-defect doubling steps from $n = 36$ to $n = 2304$. |
| [`thermal-mpo-gap-2plus1d.tex`](thermal-mpo-gap-2plus1d.tex) | [PDF](thermal-mpo-gap-2plus1d.pdf) | *The Thermal-Purity Recursion in 2+1 Dimensions*. Develops the partition-function-excess recursion for frustration-free Hamiltonians in two spatial dimensions, tested numerically on the hexagonal spin-3/2 AKLT model in [`../models/aklt/`](../models/aklt/). |

## Build

Building the notes requires a LaTeX installation with REVTeX 4.2, TikZ, and `latexmk`. From the
repository root:

```sh
mkdir -p work/latex-build
latexmk -pdf -interaction=nonstopmode -halt-on-error -file-line-error \
  -outdir=work/latex-build notes/thermal-mpo-gap-note.tex
cp work/latex-build/thermal-mpo-gap-note.pdf notes/thermal-mpo-gap-note.pdf
```

The root `.latexmkrc` adds `notes/` to `TEXINPUTS` so that `latexmk` finds `tnnotes.sty`. The `work/`
directory is ignored by git.

## Shared notation and diagram macros

The notes load [`tnnotes.sty`](tnnotes.sty), which defines shared notation macros, colors, TikZ styles,
and four general diagram building blocks so that tensor-network equations can be written inline rather
than as one-off macros:

- `\tnchain[keys]{items}`: draws a horizontal tensor network along a virtual wire, spaced automatically
  by box widths. Items are math expressions (`\cdots` produces an ellipsis). Keys: `node=`, `last=`,
  `item=`, `legs`, `ring`, `capl=`, `capr=`, `caps`, `nostubl`, `nostubr`, `name=`. Example:
  `\tnchain[node=mpo,legs,ring]{W,W,\cdots,W,W}`.
- `\tncolumn[keys]{items}`: draws a vertical tensor network along imaginary time (physical double wire).
- `prooflanes` with `\lanebound`, `\lanemove`, `\laneuse`, `\lanetag`: formats a multi-step proof
  top-to-bottom in two parallel lanes, with one box per inequality and one numbered circle per deductive
  step (a thin gray arrow into a circle marks a secondary input).
- `\sparkline{min}{max}{data}`, `\specbars{data}`, and inline word-sized glyphs (`\levels`,
  `\decayglyph`, `\fcurve`, `\overlapglyph`, `\fieldsame`, `\fieldopp`, `\stackbar`, `\torusseq`,
  `\Zt`, `\Zbox`, `\virtleg`, `\physleg`).

Only three genuinely two-dimensional diagrams are defined as standalone macros in the note's preamble
(`\Wlegs`, `\folded`, `\stackedmpo`). Use the canonical notation table below across all notes and
model documentation:

| Object | Symbol | Standard term |
| --- | --- | --- |
| Ring Hamiltonian on $L$ unit cells | $H_L = \sum_{j=1}^L h_{j,j+1} + \sum_{j=1}^L o_j$ | bond term $h$, real symmetric on-site term $o$ |
| Bond reflection-positive factorization | $-h = \sum_\alpha C_\alpha \otimes F \overline{C_\alpha} F$ | bond factors $C_\alpha$, internal cell reflection $F$ ($F = I$ for 1-site cells) |
| Spectral gap and doubling mismatch | $\Delta_L = E_1(L) - E_0(L)$, $\delta(n) = E_0(2n) - 2E_0(n)$ | spectral gap, ground-energy doubling mismatch |
| Local thermal-MPO tensor | $W_\beta^{ss'}(a,b)$ | MPO tensor ($s,s'$ physical legs, $a,b$ virtual legs) |
| Spatial trace of $W_\beta$ | $X_\beta = \sum_s W_\beta^{ss}$ | spatial transfer operator (one copy is a *transfer vertex*) |
| Eigenvalues of $X_\beta$ | $\lambda_i(\beta)$, ordered by modulus; $\lambda_j := \lambda_0(\beta_j)$ | transfer eigenvalues |
| Normalized spatial transfer power | $\sigma_n(\beta) = X_\beta^n / Z_n(\beta)$ (even $n$), weights $p_0 \ge p_1 \ge \cdots$ | normalized spatial state |
| Gibbs operator | $e^{-\beta H_n}$ (per-site shift $\kappa$ absorbed into $H_n$) | Gibbs operator |
| Periodic and open partition functions | $Z_n(\beta)$; $Y_n(\beta)$, $Y_n^{\mathrm{s}}(\beta)$ | partition function |
| Gibbs state and Boltzmann weights | $\rho_n(\beta) = e^{-\beta H_n}/Z_n(\beta)$; $w_0 \ge w_1 \ge \cdots$ | Gibbs state |
| $Z_n(2\beta)/Z_n(\beta)^2 = \operatorname{Tr}\rho_n(\beta)^2$ | $P_\tau(n,\beta)$; open chains $P_\tau^{\mathrm{o}}$, $P_\tau^{\mathrm{s}}$ | thermal purity (never "Gibbs purity" or "temporal purity") |
| $Z_{2n}(\beta)/Z_n(\beta)^2 = \operatorname{Tr}\sigma_n(\beta)^2$ | $P_{\mathrm{s}}(n,\beta)$ | spatial purity |
| $1 - P_{\mathrm{s}}(n,\beta)$ and $1 - P_\tau(n,\beta)$ | $\varepsilon_{\mathrm{s}}, \varepsilon_\tau$ at the seed $(n_0, \beta_0)$; diagonal sequence $\varepsilon_j$ | purity defects (spatial defect and thermal defect) |
| Quadratic defect recurrence map | $f(z) = z^2 / [2(1-z)^2]$, fixed point $\varepsilon_* \approx 0.12731$ | squaring map / squaring lemma |
| Cartesian spin generators and Poisson events | $G_\alpha = iS^\alpha$ ($\alpha \in \{x,y,z\}$), $(\tau,\alpha)$, measure $\mu_\beta$ | Cartesian spin generators |
| Dyadic spacetime scales | $(n_j, \beta_j) = (2^j n_0, 2^j \beta_0)$ | dyadic scales |
| Per-site energy shift | $\kappa$ (denoted $a$ in upstream `spin1` code) | per-site energy shift |

Diagram color and wire conventions are also standardized across the notes:
- **Colors:** Blue marks thermal objects ($W_\beta$, $\rho_n$) and red marks spatial objects ($X_\beta$, $\sigma_n$, $\Pi$, $R$).
- **Wires and boundaries:** Double vertical lines denote physical legs, single horizontal lines denote virtual legs, closed loops denote traces, gray circles denote boundary vectors, and dashed lattice patches denote periodic spacetime tori.
