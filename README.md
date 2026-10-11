# Thermal-MPO gap note

A concise REVTeX 4.2 note deriving a spin-1 gap bound using tensor-network
diagrams. The argument assumes the two finite-size estimates cited in the
note; their numerical certification is not reproduced here.

- Source: [outputs/thermal-mpo-gap-note.tex](outputs/thermal-mpo-gap-note.tex)
- PDF: [outputs/thermal-mpo-gap-note.pdf](outputs/thermal-mpo-gap-note.pdf)

## Build

Requires a LaTeX installation with REVTeX 4.2, TikZ, and `latexmk`.
From the repository root:

```sh
mkdir -p work/latex-build
latexmk -pdf -interaction=nonstopmode -halt-on-error -file-line-error \
  -outdir=work/latex-build outputs/thermal-mpo-gap-note.tex
cp work/latex-build/thermal-mpo-gap-note.pdf outputs/thermal-mpo-gap-note.pdf
```

The source and final PDF are versioned. The `work/` directory holds local
references, drafts, and build products and is ignored by Git.
