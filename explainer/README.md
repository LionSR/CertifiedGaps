# Spacetime bootstrap explainer

An illustrated, undergraduate-level walk through the spacetime bootstrap (two-purity doubling
criterion) and the spin-1 Haldane-gap proof in [`notes/`](../notes/), in Edward Tufte's style. It is
an exposition only: nothing here is part of any certificate.

## Contents

| file | content |
|---|---|
| `index.html` | the explainer page: imaginary time and path sums, the spacetime torus, the two purities, squaring lemma, coupling identity, reflection positivity (with the realignment test), four-move doubling, the numerical seed, and the open-chain argument of Appendices A–C, each with an interactive figure and a worked example |
| `film.html`, `film_open.html` | animation pages for the two films (main proof; open chains), drawn frame by frame by `window.draw(t)` |
| `core.js` | the exact two-purity recurrence in log space (mirrors `julia/src/recurrence/map.jl`), transverse-field Ising free-fermion purities, and the spacetime-sheet → torus renderer |
| `chapters.js`, `open_chapters.js` | the film chapters (shared by the page figures and the films) |
| `render.mjs` | renders a film page to MP4 |
| `poster.jpg`, `poster-open.jpg` | video posters |

## View

Serve the directory and open `index.html` (MathJax and EB Garamond load from public CDNs):

```sh
cd explainer && python3 -m http.server 8000   # then open http://localhost:8000
```

The page expects `film.mp4` and `open-chains.mp4` next to it. They are not committed (about
20 MB); render them with Playwright and ffmpeg:

```sh
cd explainer && npm install playwright
node render.mjs film.html film.mp4             # ~7 min of video (open film ~4 min)
node render.mjs film_open.html open-chains.mp4
```

## What is rigorous

Nothing. All numbers on the page are floating-point illustrations. The spin-1 seed numbers
(defects at $(36, 49/4)$ and the doubling chain to $(2304, 784)$) are copied from
`models/spin1/results/scalars.json` and the notes; the transverse-field Ising lab computes free-fermion
purities in `Float64`; the open-chain energy levels are a small floating-point exact diagonalization
(open spin-1 chains, $n = 4 \dots 12$, endpoint fields $h = 0$ and $3/5$) done outside the repository.

## Language note

`AGENTS.md` restricts code to Julia, Lean 4, and LaTeX. This directory is a browser-side exposition
and is written in HTML and JavaScript because that is what runs in a browser; it holds no certificate
code and nothing in `julia/`, `models/`, or `lean/` depends on it.
