# Spacetime bootstrap explainer

An illustrated, undergraduate-level walk through the spacetime bootstrap (two-purity doubling
criterion) and the spin-1 Haldane-gap proof in [`notes/`](../notes/), in Edward Tufte's style. It is
an exposition only: nothing here is part of any certificate.

## Contents

| file | content |
|---|---|
| `index.html` | Part I, the ring proof: imaginary time and path sums, the spacetime torus, the two purities, squaring lemma, coupling identity, reflection positivity (with the realignment test), four-move doubling, and the numerical seed, each with an interactive figure and a worked example |
| `open.html` | Part II, the open-chain argument of Appendices A–C: boundary vectors, the overlap, the three moves per scale, and the odd-length gap |
| `page.css`, `page.js` | the style and runtime both pages share: drawing helpers, theme colors, TeX on canvas, and the figure players |
| `film.html`, `film_open.html` | animation pages for the two films (main proof; open chains); each only names its chapter list and title |
| `film.js` | the shared film player: TeX preloading, captions, the beat timeline, and the frame-deterministic `window.draw(t)` |
| `narration.js` | the spoken narration of both films, one line per caption beat, with symbols spelled out for speech |
| `core.js` | the exact two-purity recurrence in log space (mirrors `julia/src/recurrence/map.jl`), transverse-field Ising free-fermion purities, and the spacetime-sheet → torus renderer |
| `chapters.js`, `open_chapters.js` | the film chapters (shared by the page figures and the films); `chapters.js` also holds the drawing kit and the helpers both files use |
| `render.mjs` | renders a film page to MP4 with its narration track |
| `poster.jpg`, `poster-open.jpg` | video posters |

## View

Serve the directory and open `index.html` or `open.html` (MathJax and EB Garamond load from public CDNs):

```sh
cd explainer && python3 -m http.server 8000   # then open http://localhost:8000
```

The pages expect `film.mp4` and `open-chains.mp4` next to them. They are too large for git, so they
live as assets of the [`films` release](https://github.com/LionSR/CertifiedGaps/releases/tag/films),
and the site deploy downloads them from there. To view them locally, download them
(`gh release download films --dir explainer`) or render them with Playwright, ffmpeg, and the
[Kokoro-82M](https://huggingface.co/hexgrad/Kokoro-82M) speech model (run locally through
`kokoro-js`; its ONNX weights, about 300 MB, download on first use). Rendering both films takes about
17 minutes on a laptop once the narration clips are cached:

```sh
cd explainer && npm install --no-save playwright kokoro-js
npx playwright install chromium                  # once, if Playwright has no browser yet
node render.mjs film.html film.mp4               # narrated main film
node render.mjs film_open.html open-chains.mp4   # narrated open-chain film
node render.mjs film.html film.mp4 --mute        # silent, on the original caption clock
```

Each caption beat is read from `narration.js` in the `af_heart` voice. A beat whose line runs longer
than its caption plays up to 25% slower and then holds its last frame until the line ends. Clips are
cached in `.voice/` by a hash of their text, so editing one line re-synthesizes only that line.

After re-rendering, replace the release assets; the next deploy (any push to `main`) serves them:

```sh
gh release upload films film.mp4 open-chains.mp4 --clobber
```

## AI use

Claude Opus 5.5 (Anthropic), run in Claude Code, wrote the text, figures, and code of both pages
and the films, from prompts by Sirui Lu. The Kokoro-82M speech model reads the narration. The proofs
are from the OpenAI preprints. Each page states this below its title.

## What is rigorous

Nothing. All numbers on the page are floating-point illustrations. The spin-1 seed numbers
(defects at $(36, 49/4)$ and the doubling chain to $`(2304, 784)`$) are copied from
`models/spin1/results/scalars.json` and the notes; the transverse-field Ising lab computes free-fermion
purities in `Float64`; the open-chain energy levels are a small floating-point exact diagonalization
(open spin-1 chains, $n = 4 \dots 12$, endpoint fields $h = 0$ and $3/5$) done outside the repository.
