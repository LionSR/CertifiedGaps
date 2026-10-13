# Notes

REVTeX 4.2 notes. Each `.tex` has its compiled `.pdf` beside it.

| note | PDF | content |
|---|---|---|
| [`thermal-mpo-gap-note.tex`](thermal-mpo-gap-note.tex) | [PDF](thermal-mpo-gap-note.pdf) | *A Diagrammatic Proof of the Haldane Gap via Spacetime Tensor Networks*. Main text: a uniform gap Δ_L > (log 20)/784 ≈ 0.003821 for the periodic spin-1 antiferromagnetic Heisenberg chain, from an exact continuous-time thermal MPO and the spacetime bootstrap, a scale-doubling iteration on two purities. Appendices A–C: open chains with equal endpoint fields h = 3/5, Δ_n > (log 10)/392 for odd n ≥ 1921. |
| [`thermal-mpo-gap-scope.tex`](thermal-mpo-gap-scope.tex) | [PDF](thermal-mpo-gap-scope.pdf) | *Scope of the Spacetime Bootstrap* (appendices D–H). Which Hamiltonians the method covers: bond reflection positivity makes the spatial transfer operator self-adjoint, and bounds on the two purities at one scale then give a uniform gap. |
| [`thermal-mpo-gap-certificate.tex`](thermal-mpo-gap-certificate.tex) | [PDF](thermal-mpo-gap-certificate.pdf) | *How the Initial Estimate Is Certified* (appendices I–M). In tensor-network diagrams: twisted partition functions, symmetry-sector projectors, polynomial filters, one matrix-product trial vector, and six unequal-defect doublings that certify the estimate at (2304, 784). |
| [`thermal-mpo-gap-2plus1d.tex`](thermal-mpo-gap-2plus1d.tex) | [PDF](thermal-mpo-gap-2plus1d.pdf) | *The thermal-purity recursion in 2+1 dimensions*: the partition-function-excess recursion for frustration-free Hamiltonians; tested numerically on the hexagonal AKLT model in [`../models/aklt/`](../models/aklt/). |

## Build

Requires a LaTeX installation with REVTeX 4.2, TikZ, and `latexmk`.
From the repository root:

```sh
mkdir -p work/latex-build
latexmk -pdf -interaction=nonstopmode -halt-on-error -file-line-error \
  -outdir=work/latex-build outputs/thermal-mpo-gap-note.tex
cp work/latex-build/thermal-mpo-gap-note.pdf outputs/thermal-mpo-gap-note.pdf
```

The root `.latexmkrc` adds `outputs/` to `TEXINPUTS`, so the build finds `tnnotes.sty`.
`work/` is not version controlled.

## Shared notation

The note loads `outputs/tnnotes.sty`, which holds the shared notation macros,
colours and TikZ styles, and four general building blocks; diagrams are written
inline with them rather than as one-off macros:

- `\tnchain[keys]{items}`: a horizontal network on a virtual wire, spaced by the
  real box widths. Items are math; `\cdots` is an ellipsis. Keys: `node=`,
  `last=`, `item=`, `legs`, `ring`, `capl=`, `capr=`, `caps`, `nostubl`,
  `nostubr`, `name=`. Example: `\tnchain[node=mpo,legs,ring]{W,W,\cdots,W,W}`.
- `\tncolumn[keys]{items}`: the same along imaginary time (physical double wire).
- `prooflanes` with `\lanebound`, `\lanemove`, `\laneuse`, `\lanetag`: a proof
  read top to bottom in two lanes, one box per bound with its full arguments and
  one numbered circle per move; a thin grey arrow into a circle is a second input.
- `\sparkline{min}{max}{data}`, `\specbars{data}` and the word-sized glyphs
  (`\levels`, `\decayglyph`, `\fcurve`, `\overlapglyph`, `\fieldsame`,
  `\fieldopp`, `\stackbar`, `\torusseq`, `\Zt`, `\Zbox`, `\virtleg`, `\physleg`).

Only three genuinely two-dimensional diagrams remain as macros in the note's
preamble (`\Wlegs`, `\folded`, `\stackedmpo`). Follow this table when editing.

| Object | Symbol | Term |
| --- | --- | --- |
| Local thermal-MPO tensor | $W_\beta^{ss'}(a,b)$ | MPO tensor; $s,s'$ physical, $a,b$ virtual legs |
| $\sum_s W_\beta^{ss}$ | $X_\beta$ (as in both source papers) | spatial transfer operator; one copy is a *transfer vertex* |
| Eigenvalues of $X_\beta$ | $\lambda_i(\beta)$, ordered by modulus; $\lambda_j:=\lambda_0(\beta_j)$ | — |
| Eigenvalues of $\sigma_n$ / any normalized spectrum | $p_0\ge p_1\ge\cdots$ | — |
| Gibbs operator | $e^{-\beta H_n}$ (no separate symbol; the shift $\kappa$ is absorbed into $H$) | — |
| Periodic / open partition functions | $Z_n(\beta)$; $Y_n(\beta)$, $Y_n^{\mathrm s}(\beta)$ | — |
| Gibbs state and its weights | $\rho_n(\beta)$; $w_0\ge w_1\ge\cdots$ | — |
| $Z_n(2\beta)/Z_n(\beta)^2=\operatorname{Tr}\rho_n^2$ | $P_\tau(n,\beta)$; open chains $P_\tau^{\mathrm o}$, $P_\tau^{\mathrm s}$ | thermal purity (not "Gibbs purity") |
| $Z_{2n}(\beta)/Z_n(\beta)^2$ | $P_{\mathrm s}(n,\beta)$ | spatial purity |
| $1-P$ | — | purity defect |
| Periodic purity-defect sequence | $\varepsilon_j$ | — |
| Squaring bound | $f(z)=z^2/[2(1-z)^2]$ | squaring lemma |
| Spin generators, bond events | $G_\alpha=iS^\alpha$, $(\tau,\alpha)$, measure $\mu_\beta$ | — |
| Scales | $\beta_j$, $n_j$ | — |
| Energy shift per site (open-chain appendix) | $\kappa$ (called $a$ in the paper) | — |

Diagram conventions are shared too. Blue marks thermal objects ($W$, $\rho$) and red marks spatial ones ($X$, $\sigma$, $\Pi$, $R$).
Double vertical lines are physical legs, horizontal lines are virtual legs, and closed loops are traces.
Grey circles are boundary vectors. Dashed lattice patches are periodic spacetime tori.
