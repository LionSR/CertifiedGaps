# LaTeX Notes

This directory contains the REVTeX 4.2 notes that develop the spacetime bootstrap, its scope across
reflection-positive spin and ladder models, and the diagrammatic certification of the finite seed
inputs. Each `.tex` source file has its compiled `.pdf` alongside it. The companion notes
`thermal-mpo-gap-scope.tex` (Appendices D–H) and `thermal-mpo-gap-certificate.tex` (Appendices I–M)
continue the appendix numbering of `thermal-mpo-gap-note.tex` (Appendices A–C).

| note | PDF | content |
|---|---|---|
| [`thermal-mpo-gap-note.tex`](thermal-mpo-gap-note.tex) | [PDF](thermal-mpo-gap-note.pdf) | *A Diagrammatic Proof of the Haldane Gap via Spacetime Tensor Networks*. Main text: proves a uniform spectral gap $\Delta_L > (\log 20)/784 \approx 0.003821$ for the periodic spin-1 antiferromagnetic Heisenberg chain for all even $L \ge 2304$, using an exact continuous-time thermal MPO and a scale-doubling iteration on spatial and thermal purities. Appendices A–C treat open chains with equal endpoint fields $h = 3/5$, proving $\Delta_n > (\log 10)/392$ for odd $n \ge 1921$. |
| [`spin1-gap-tn-diagrams.tex`](spin1-gap-tn-diagrams.tex) | [PDF](spin1-gap-tn-diagrams.pdf) | *A Short Diagrammatic Explanation of the Spin-1 Haldane Gap Proof* (S. Lu and GPT-6-Astra, October 7, 2026). Self-contained two-page version of the main-text argument: from the two finite-size estimates of the OpenAI periodic preprint, derives $\Delta_L > (\log 20)/784$ for every even $L \ge 2304$. Standalone source (no `tnnotes.sty`). |
| [`thermal-mpo-gap-scope.tex`](thermal-mpo-gap-scope.tex) | [PDF](thermal-mpo-gap-scope.pdf) | *Scope of the Spacetime Bootstrap* (Appendices D–H). Formalizes the general criterion for translation-invariant chains and extends it to multi-site unit cells with internal reflection $F$ (`RPRing`), positive cell-transfer factorizations $X_\beta = KK^* \succeq 0$ for odd seed lengths, the asymmetric $(\varepsilon_{\mathrm{s}}, \varepsilon_\tau)$ recurrence basin, and half-length seeds; classifies bilinear, biquadratic, staggered-field, anisotropic, dimerized, and two-leg ladder Hamiltonians; and summarizes the computer-verified gap certificates across eight model families (`CertifiedGaps.VerifiedSpectra`, coupling-interval perturbations, all-$`\beta_0`$ exclusions, exact rational chain certificates, and the unconditional Lean proof of the longitudinal Ising seed at $`(g,h) = (1/4, 1/20)`$). |
| [`thermal-mpo-gap-certificate.tex`](thermal-mpo-gap-certificate.tex) | [PDF](thermal-mpo-gap-certificate.pdf) | *How the Initial Estimate Is Certified* (Appendices I–M). Presents the finite seed certificate at $(n_0, \beta_0) = (2304, 784)$ in tensor-network diagrams: twisted partition functions, symmetry-sector projectors, degree-12 polynomial filters, a matrix-product trial vector at $n = 72$, and six unequal-defect doubling steps from $(36, 49/4)$ (including temperature lifting on the $b_0 = 49/4$ sector moments giving $\Delta_L > 0.0729$ for even $L \ge 36$ and the untwisted one-trace specialization), followed by smaller seed inputs and half-chain comparison operators for formal verification (`models/lean_seeds`), $\mathrm{SU}(2)$ reduced-block trial certificates (`models/spin1`, `models/spin2`), and certified thermal MPO compression (`CertifiedGaps.CertifiedMPO`: Taylor and Araki–Lieb–Thirring majorants, rational matrix gauges, and three-layer trace corrections with backward error propagation). |
| [`thermal-mpo-gap-2plus1d.tex`](thermal-mpo-gap-2plus1d.tex) | [PDF](thermal-mpo-gap-2plus1d.pdf) | *The Thermal-Purity Recursion in 2+1 Dimensions*. Develops the partition-function-excess recursion for frustration-free Hamiltonians in two spatial dimensions, tested numerically on the hexagonal spin-3/2 AKLT model in [`../models/aklt/`](../models/aklt/). |

## References

[`references.bib`](references.bib) is the concise reference list for the arXiv version of the paper
(26 BibTeX entries). Every entry was checked against its Crossref DOI record or arXiv page on
2026-10-16. The notes above do not use it yet; they keep their own bibliographies.
[`references.html`](references.html) shows the same list with a small figure per entry of what that
paper did, a line on why to cite it, a suggested disclaimer, and the papers left out; open it in a
browser.

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
  `\decayglyph`, `\fcurve`, `\overlapglyph`, `\fieldsame`, `\fieldopp`, `\fieldstag`, `\stackbar`,
  `\torusseq`, `\Zt`, `\Zbox`, `\virtleg`, `\physleg`).

Only three genuinely two-dimensional diagrams are defined as standalone macros in the note's preamble
(`\Wlegs`, `\folded`, `\stackedmpo`). Use the canonical notation table below across all notes and
model documentation:

| Object | Symbol | Standard term |
| --- | --- | --- |
| Ring Hamiltonian on $L$ unit cells | $H_L = \sum_{j=1}^L h_{j,j+1} + \sum_{j=1}^L o_j$ | bond term $h$, real symmetric on-site term $o$ |
| Bond reflection-positive factorization | $-h = \sum_\alpha C_\alpha \otimes F \overline{C_\alpha} F$ | bond factors $C_\alpha$, internal cell reflection $F$ ($F = I$ for 1-site cells) |
| Spectral gap | $\Delta_L = E_1(L) - E_0(L)$ | spectral gap |
| Local thermal-MPO tensor | $W_\beta^{ss'}(a,b)$ | MPO tensor ($s,s'$ physical legs, $a,b$ virtual legs) |
| Spatial trace of $W_\beta$ | $X_\beta = \sum_s W_\beta^{ss}$ | spatial transfer operator (one copy is a *transfer vertex*) |
| Eigenvalues of $X_\beta$ | $\lambda_i(\beta)$, ordered by modulus; $\lambda_j := \lambda_0(\beta_j)$ | transfer eigenvalues |
| Normalized spatial transfer power | $\sigma_n(\beta) = X_\beta^n / Z_n(\beta)$ (even $n$), weights $p_0 \ge p_1 \ge \cdots$ | normalized spatial state |
| Gibbs operator | $e^{-\beta H_n}$ (per-site shift $\kappa$ absorbed into $H_n$) | Gibbs operator |
| Periodic and open partition functions | $Z_n(\beta)$; $Y_n(\beta)$, $Y_n^{\mathrm{s}}(\beta)$ | partition function |
| Gibbs state and Boltzmann weights | $\rho_n(\beta) = e^{-\beta H_n}/Z_n(\beta)$; $w_0 \ge w_1 \ge \cdots$ | Gibbs state |
| $`Z_n(2\beta)/Z_n(\beta)^2 = \mathrm{Tr}\,\rho_n(\beta)^2`$ | $P_\tau(n,\beta)$; open chains $P_\tau^{\mathrm{o}}$, $P_\tau^{\mathrm{s}}$ | thermal purity (never "Gibbs purity" or "temporal purity") |
| $`Z_{2n}(\beta)/Z_n(\beta)^2 = \mathrm{Tr}\,\sigma_n(\beta)^2`$ | $P_{\mathrm{s}}(n,\beta)$ | spatial purity |
| $1 - P_{\mathrm{s}}(n,\beta)$ and $1 - P_\tau(n,\beta)$ | $\varepsilon_{\mathrm{s}}, \varepsilon_\tau$ at the seed $(n_0, \beta_0)$; diagonal sequence $\varepsilon_j$ | purity defects (spatial defect and thermal defect) |
| Quadratic defect recurrence map | $f(z) = z^2 / [2(1-z)^2]$, fixed point $\varepsilon_* \approx 0.12731$ | squaring map / squaring lemma |
| Cartesian spin generators and Poisson events | $G_\alpha = iS^\alpha$ ($`\alpha \in \{x,y,z\}`$), $(\tau,\alpha)$, measure $\mu_\beta$ | Cartesian spin generators |
| Dyadic spacetime scales | $(n_j, \beta_j) = (2^j n_0, 2^j \beta_0)$ | dyadic scales |
| Per-site energy shift | $\kappa$ (denoted $a$ in upstream `spin1` code) | per-site energy shift |

Diagram color and wire conventions are also standardized across the notes:
- **Colors:** Blue marks thermal objects ($W_\beta$, $\rho_n$) and red marks spatial objects ($X_\beta$, $\sigma_n$, $\Pi$, $R$).
- **Wires and boundaries:** Double vertical lines denote physical legs, single horizontal lines denote virtual legs, closed loops denote traces, gray circles denote boundary vectors, and dashed lattice patches denote periodic spacetime tori.
