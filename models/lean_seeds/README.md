# Toward a Lean-checked seed

The Lean 4 gap proofs in [`lean/`](../../lean) formalize the reflection-positivity transfer
operator and the all-scale doubling criterion, leaving each model's finite seed inequality as an
explicit `sorry` hypothesis in the model file. This directory contains five studies investigating
how to discharge those finite seed inequalities inside the Lean kernel using `decide +kernel`
alone (with no `native_decide` and no additional axioms): four studies targeting the spin-1
Heisenberg chain and one end-to-end pipeline on the computationally cheapest certified model,
the longitudinal-field Ising chain. Throughout these studies, Julia serves as an untrusted
witness generator, while the verified Lean declarations live in the opt-in library `SeedKernel`
of the repository's Lake package (`lake build SeedKernel`).

**One unconditional gap.** The Ising study closes: `SeedKernel.Ising.ising_gap` proves, with no
hypothesis and no `sorry`, that the Ising ring at $g = 1/4$, $h = 1/20$ has spectral gap
$\Delta_L > \log(13480)/(97/4) \approx 0.3921$ for every even $L \ge 6$. The spin-1 studies remain
exact but conditional on open bridge lemmas or spectral counting anchors:

| Study | Core question | Best result and verification status |
|---|---|---|
| `seed_design` | Which finite inputs does the spin-1 seed actually require? | Design C: $\Delta_L > 0.039036761$ for all even $L \ge 60$ from even $n \le 12$, real twists ($I, P$), and trace error $10^{-3}$ (exact rational replay of `ChainCertificate.Valid`; conditional on trace tables and analytic lemmas) |
| `spectrum` | How many low-lying states must be resolved, to what precision, and at what kernel cost? | $26$ ($I$) and $34$ ($P$) states below the tail cutoff at $n = 12, \beta = 49/4$; degree-16 Q48 Chebyshev filtering meets all residual targets (`Float64` cost study, non-certifying) |
| `anchor` | Can a local comparison operator anchor the low-energy eigenvalue count? | Exact local PSD comparison via rational $\mathrm{LDL}^\top$; reduces $S^z = 0$ comparison count at $n = 12$ from $191 \to 121 \to 49$ (binned tail) against $8$ physical levels (exact counts; one-shot anchor does not close without continuation) |
| `spin1_kernel` | Does base-3 packed integer arithmetic fit within the Lean kernel? | Yes: $n = 12$ dense matvec + dot checks in $20\text{--}26\text{ s}$; power-filter Ritz and residual witnesses for $H_6, H_8, H_{10}$ kernel-checked (arithmetic inequalities only) |
| `ising` | Can we run the complete pipeline on the cheapest certified model? | **Proved**: seed of `IsingLongitudinal` at $g = 1/4, h = 1/20, (n_0, \beta_0) = (6, 97/4)$ and $\Delta_L > \log(13480)/(97/4) \approx 0.3921$ for every even $L \ge 6$ (`ising_gap`; standard axioms only) |

## Model

- **Spin-1 Heisenberg chain** (`seed_design`, `spectrum`, `anchor`, `spin1_kernel`):
  ```math
  H_n = \sum_{j=1}^n \mathbf{S}_j \cdot \mathbf{S}_{j+1}
  ```
  on a periodic ring of $n$ sites (lengths count sites), shifted by $\kappa = 700741/500000$ per
  bond (denoted $a$ in the upstream source code), with the closing bond twisted by a rotation
  about the $z$-axis by angle $0$ ($I$), $\pi$ ($P$), or $2\pi/3$ ($C$). The numerical and
  Lean-kernel studies work in the $S^z$ product basis indexed by $\sum_{j=0}^{n-1} (m_j + 1)3^j$,
  where the two-site bond operator $h = S^z \otimes S^z + \frac{1}{2}(S^+ \otimes S^- + S^- \otimes S^+)$
  has integer entries and eigenvalues $-2, -1, 1$. Source 1 evaluates $\beta = 49/4$ on lengths
  $n = 4\text{--}12$ with twists $I, P$; source 2 evaluates $\beta = 21/2$ on lengths
  $n = 6, 8, 10$ with twists $I, P, C$ (`models/spin1/results/source{1,2}.json`).
- **Transverse- and longitudinal-field Ising chain** (`ising`):
  ```math
  H_L = -\sum_{j=1}^L Z_j Z_{j+1} - g \sum_{j=1}^L X_j - h \sum_{j=1}^L Z_j
  ```
  at $g = 1/4, h = 1/20$, matching `CertifiedGaps.IsingLongitudinal.model (1/4) (1/20)`. The
  kernel checkers scale the Hamiltonian to the integer matrix $A = 20 H_L$ in the computational
  ($\sigma^z$) basis; `test/ising.jl` verifies that its periodic-bond convention at $L = 2$
  matches `RPRing.hamiltonian`.

## Reflection positivity

These studies focus on the finite seed inequalities rather than reflection positivity, which is
already formalized in [`lean/`](../../lean):
- For the spin-1 Heisenberg chain, reflection positivity is proved in the Cartesian basis in
  `lean/CertifiedGaps/Heisenberg/Main.lean` ($-h = \sum_\alpha C_\alpha \otimes \overline{C_\alpha}$
  with $C_\alpha = i S^\alpha$); the unitary change-of-basis bridge from that Cartesian matrix to
  the $S^z$-basis integer matrices used in `SeedKernel` remains an open formalization task.
- For the longitudinal-field Ising chain, reflection positivity holds in the standard Pauli-$`Z`$
  basis with single bond factor $C_1 = Z$ ($-h_{\mathrm{bond}} = Z \otimes Z$) and real diagonal
  on-site term $o = -g X - h Z$, checked directly in `lean/CertifiedGaps/Models/IsingLongitudinal.lean`.

## Result

### 1. `seed_design`: smaller finite seed inputs for spin-1

Three simplified seed designs were found and verified in exact rational arithmetic down to
`ChainCertificate.Valid` (`results/seed_design/certificate.json`, using 9-step doubling chains on
the rational rounding grid $10^{-512}$). These bounds are **conditional** on certifying the
short-ring trace tables to the stated absolute error $\epsilon$ and formalizing the analytic
sector, polynomial-filter, and temperature-lifting lemmas ([`METHOD.md`](METHOD.md) §1). Note
that in `seed_design`, the spatial and thermal seed defects are denoted $(p, q) = (\varepsilon_{\mathrm{s}}, \varepsilon_\tau)$.

| Design | Thermal trace inputs (lengths; twists; abs. error $\epsilon$) | $\beta_{\mathrm{table}} \to \beta_0$ | Moment $k$; trial $N$ | $\varepsilon_{\mathrm{s}}$ ($p$) | $\varepsilon_\tau$ ($q$) | Certified gap $\Delta_L >$ | Lengths covered |
|---|---|---|---|---:|---:|---:|---|
| A (`small_simple`) | $`n \in \{6, 8, 10\}`$; $I, P, C$; $3/10000$ | $21/2 \to 105/4$ | $30$; $120$ | $0.124048525$ | $0.125580476$ | $0.005966152$ | even $L \ge 60$ |
| B (`small_relaxed`) | $`n \in \{6, 8, 10\}`$; $I, P, C$; $4/10000$ | $21/2 \to 525/16$ | $30$; $144$ | $0.117346685$ | $0.128925715$ | $0.006309519$ | even $L \ge 72$ |
| C (`real_even`) | $`n \in \{4, 6, 8, 10, 12\}`$; $I, P$; $1/1000$ | $49/4 \to 49/2$ | $30$; $120$ | $0.047914768$ | $0.112578548$ | $0.039036761$ | even $L \ge 60$ |

All three designs improve on the published source-1 bound $\log(20)/784 \approx 0.00382$ while
eliminating odd chain lengths and loosening the required trace tolerance from $3 \times 10^{-5}$
to $3 \times 10^{-4}\text{--}10^{-3}$. Design B's thermal defect $\varepsilon_\tau$ slightly
exceeds the diagonal threshold $\varepsilon_* \approx 0.12731$, but its unequal-defect doubling
chain converges into the admissible region. The exploratory scan (`diagnostics.json`) also
established three sharp negative results within the tested filter family: restricting to lengths
$n \le 8$ fails even at zero trace error (best $\max(\varepsilon_{\mathrm{s}}, \varepsilon_\tau) \approx 0.5288$
across $132{,}075$ temperature-lifting configurations); restricting to $n \le 10$ without the
complex twist $C$ fails ($\approx 0.2956$); and $n \le 10$ with twist $C$ at error $10^{-3}$
narrowly misses the admissible threshold ($\approx 0.1338$).

### 2. `spectrum`: low-energy state counts, filter precision, and kernel cost models

| Source | $\beta_{\mathrm{table}}$ | Rings included | States below tail cutoff at $n_{\max}$ ($I$ / $P$ / $C$) | Adaptive Q48 $Hv$ products | Projected Lean vector + Gram time |
|---|---:|---|---|---:|---|
| Source 1 | $49/4$ | $n = 4\text{--}12$, $I/P$ | $26$ / $34$ / — at $n = 12$ | $34{,}388$ | $38\text{--}83\text{ min}$ |
| Source 2 | $21/2$ | $n = 6, 8, 10$, $I/P/C$ | $22$ / $30$ / $26$ at $n = 10$ | $17{,}800$ | $5\text{--}11\text{ min}$ |

Key findings from the floating-point spectral and fixed-point filtering experiments:
- **Chebyshev filtering precision:** 48-bit fixed-point (Q48) Chebyshev filtering with degree-16
  cycles and periodic Rayleigh–Ritz orthogonalization meets every sector's eigenvalue residual
  target. Increasing the polynomial degree to 32 in Q48 without orthogonalization causes
  stagnation at $2.99 \times 10^{-6}$ for $(n = 12, I, S^z = 0)$, whereas degree-16 Q48 reaches
  $3.29 \times 10^{-8}$ and Q112 reaches $2.34 \times 10^{-9}$.
- **Basis layout:** Directly packing compressed $S^z$ sectors destroys fixed shift distances for
  nearest-neighbor spin flips, requiring $3{,}816$ distinct displacement masks at $(n = 12, S^z = 0)$
  and $13{,}620$ at $S^z = 2$ instead of $2n = 24$ shifts in the full $3^n$ ternary basis. Full
  ternary packing is therefore substantially faster in the Lean kernel.
- **Naive anchor gap:** Cutting the ring into two open half-chains with the bare bound $h \ge -2I$
  leaves $242$ comparison levels in the $S^z = 0$ sector below the $n = 12$ tail cutoff $\rho$,
  compared with only $10$ physical states on the ring.

### 3. `anchor`: exact local comparison operators for spectral counting

To tighten the two-half-chain comparison $H \ge H_0 = A_l \otimes I + I \otimes A_l + 2cI$
(where $l = n/2$), we add $\mathrm{SU}(2)$-invariant boundary corrections
$A_l = H_{\mathrm{open}, l} + B_{1,2} + B_{l-1,l}$ with $B = u h + v h^2$. The operator inequality
$H - H_0 \ge 0$ holds whenever the four-site cut remainder
```math
R = h_{23} - B_{12} - B_{34} - cI
```
is positive semidefinite. This $81 \times 81$ condition is verified in exact rational arithmetic
by block $\mathrm{LDL}^\top$ factorization in every $S^z$ sector (`results/anchor/certificate.json`),
and the half-chain spectra of $A_l$ (sector dimensions at most $141$ at $l = 6$) are enclosed
rigorously to radius $\le 2.8 \times 10^{-7}$.

| $n$ | Count-optimized $(u, v, c)$ | $S^z = 0$ count below uniform cutoff: bare / ground-opt / count-opt | Binned-tail $S^z = 0$ count (all-$`S^z`$ total) | Physical $S^z = 0$ states |
|---:|---|---|---|---:|
| $8$ | $`(-3/10,\; -1/10,\; -5076841/2500000)`$ | $44$ / $63$ / $34$ | $28$ ($98$) | $2$ |
| $10$ | $`(-1/5,\; 0,\; -5538721/2500000)`$ | $86$ / $115$ / $62$ | $38$ ($140$) | $5$ |
| $12$ | $`(-1/5,\; 1/20,\; -48533/20000)`$ | $191$ / $212$ / $121$ | $49$ ($183$) | $8$ |

Replacing the coarse dimension tail $(3^n - k)e^{-\beta \rho}$ by a 32-bin histogram of the exact
comparison spectrum $H_0$ lowers the initial $S^z = 0$ state count at $n = 12$ from $121$ to
$49$. However, because the comparison count $K$ still exceeds the number of low-energy physical
states $k$, a one-shot comparison count does not immediately yield ordered lower bounds on the
first $k$ eigenvalues (as illustrated by the elementary counterexample $H = \mathrm{diag}(-2, -1, 10)$
with trial subspace $\mathrm{span}(e_2, e_3)$ and $`\rho = 0 \le \lambda_3(H)`$). Bridging
from $K$ comparison states to $k$ physical states via homotopy continuation
$H(t) = H_0 + t(H - H_0)$ succeeded numerically across all $26$ active $(n, g, S^z)$ sectors in
$658$ total stages ($\approx 40\text{ min}$ in Julia), but a naive Lean kernel translation would
require an estimated $80\text{--}97$ kernel hours at $n = 12$.

### 4. `spin1_kernel`: base-3 packed `Nat` arithmetic in the Lean kernel

The modules `lean/SeedKernel/Spin1/H{6,8,10}.lean` (using `Std` only, with no Mathlib imports)
verify integer Rayleigh–Ritz and residual inequalities (`ritzCheck`) inside the Lean kernel via
`decide +kernel` on the moments $G = v^\top v$, $A = v^\top F v$, and $C = (Fv)^\top (Fv)$ of the
Marshall-transformed power-filtered trial state $`v = (nI - H_{\mathrm{M}})^q e_{\text{N\'eel}}`$:

| $n$ | Filter degree $q$ | Lane width (bits) | Kernel-checked Ritz witness $\theta \le$ | Kernel-checked residual witness $\le$ |
|---|---:|---:|---:|---:|
| $6$ | $24$ | $262$ | $-8495447 / 10^6$ | $274967 / 10^6$ |
| $8$ | $32$ | $345$ | $-11179504 / 10^6$ | $279904 / 10^6$ |
| $10$ | $16$ | $222$ | $-13548710 / 10^6$ | $632355 / 10^6$ |

On an Apple M5 machine (`results/spin1_kernel/benchmarks.tsv`), a single 128-bit-lane dense
matrix-vector product plus dot product checks in $0.5\text{ s}$ ($n = 6, 8$), $1.3\text{ s}$
($n = 10$), and $26\text{ s}$ ($n = 12$) per fresh Lean process, or $9.6\text{ s}$ at $n = 12$
when using hoisted masks and division-free bit-reversal (`Fast12`). A full audit of the pinned
Mathlib commit (`c55e6e78`) confirms that while Hermitian spectral diagonalization and Jensen's
inequality are available, Courant–Fischer min–max, Weyl's inequalities, and Cauchy interlacing
must be formalized from scratch.

### 5. `ising`: end-to-end prototype on the longitudinal Ising chain

For the longitudinal-field Ising ring at $g = 1/4, h = 1/20$ with seed $(n_0, \beta_0) = (6, 97/4)$,
only the two lowest eigenstates (the near-degenerate ferromagnetic doublet) lie below the
classical domain-wall threshold $\rho_n = (82 - 26n)/20$. In `lean/SeedKernel/Ising/`, the Lean
kernel checks—with scale $S = 2^{52}$ and generated integer eigenvectors $v_1, v_2$ and candidate
eigenvalues $\lambda_1, \lambda_2$ for $A = 20H_n$ at $n = 6$ and $n = 12$—that:
1. $|v_i^\top v_i - S^2| \le 10^{-12} S^2$ and $|v_1^\top v_2| \le 10^{-12} S^2$,
2. $`\sum_{i=1}^2 \|S A v_i - \lambda_i v_i\|^2 \le 10^{-18} S^4`$,
3. $\lambda_1 + S < \lambda_2$ and $\lambda_2 + S < (82 - 26n) S$, and
4. exactly two classical spin configurations have diagonal energy below $(82 - 26n)/20 + n g$.

Both packed-lane (`Packed.lean`) and direct chunked (`GeneratedChunks.lean`, split into 64 blocks
at $n = 12$ with a machine-checked range-coverage theorem) verifiers are implemented, along with
negative tests confirming rejection of corrupted witnesses.

| Ring size $n$ | Hilbert-space dimension $2^n$ | Candidate $E_0$ | Candidate $E_1$ | Classical anchor $\rho_n$ | Normalized residual |
|---|---:|---:|---:|---:|---:|
| $6$ sites | $64$ | $-6.391774806847206$ | $-5.796599647006008$ | $-3.7$ | $5.97 \times 10^{-15}$ |
| $12$ sites | $4096$ | $-12.78354960854822$ | $-11.59319932209873$ | $-11.5$ | $1.49 \times 10^{-13}$ |

`lean/SeedKernel/Ising/Seed.lean` turns these checks into `first_point_seed`: exactly the first of
the thirty conjunctions of `IsingLongitudinal.seeds`
($`1 - 9064 \times 10^{-9} \le \mathrm{spatialPurity}\; 6\; (97/4)`$ and
$`1 - 3141 \times 10^{-8} \le \mathrm{thermalPurity}\; 12\; (97/4)`$), with no hypothesis, and then
`ising_gap`. The steps are:
- a residual bound gives an eigenvalue within $10^{-10}$ of each candidate (`Ising/Spectral`);
- $`\mathrm{Re}\, \langle x, Hx \rangle \ge \sum_\sigma (E(\sigma) - L/4) |x(\sigma)|^2`$ and the
  classical count give at most two eigenvalues below the anchor (counting by comparison);
- the checker's integer sums are identified with the matrix `IsingLongitudinal.model (1/4) (1/20)`
  (`Ising/Hamiltonian`, `Ising/Bridge`);
- finally
  $e^{-\beta(\theta_0 + r)} \le Z \le e^{-\beta(\theta_0 - r)}(1 + e^{-\beta(\theta_1 - \theta_0)} + 2^L e^{-\beta(\rho - \theta_0 + r)})$
  with $e > 2.7182818283$ bounds both purities.

The proof does not use the `sorry` in `IsingLongitudinal.seeds`; the earlier conditional bridge from
four $\log Z$ enclosure hypotheses is retired. A broader 43-target cost inventory
(`results/ising/costs.json`, `Float64`) across all seven certified model families shows that while
small-ring Ising seeds take seconds, naive $k \times k$ Gram checks on $65{,}536$- and
$531{,}441$-dimensional seeds at smaller $\beta_0$ would require hours to $\sim 1{,}000$ hours
without symmetry or cluster reductions.

### Open formalization obligations

Closing an unconditional kernel-checked spin-1 seed requires completing five components:
1. **Spectral counting anchor:** a kernel-checkable certificate that $\rho \le \lambda_{k+1}(H)$ for
   the spin-1 rings (via a sharper comparison operator, fewer continuation stages, or
   Lehmann–Goerisch bounds with a priori input); this is the bottleneck of the whole route.
2. **Ordered spectral linear algebra in Lean:** Lehmann bounds and ordered-eigenvalue statements
   (absent from Mathlib). The residual bound and counting by comparison are now proved
   (`Ising/Spectral.lean`) and apply to any Hermitian matrix.
3. **Packing soundness:** proofs of the specification propositions (`DotCorrect`, `HApplyCorrect`,
   `FilterCorrect`, etc.) in `lean/SeedKernel/Obligations.lean` and `Spin1/Obligations.lean`, and
   the lane decoding of `Ising/Packed.lean`; nothing assumes them.
4. **Hamiltonian basis bridges:** identification of the integer matrices checked in `SeedKernel`
   with the Lean Hamiltonians (Cartesian-to-$`S^z`$ unitary and boundary twists for spin-1); done for
   Ising (`Ising/Bridge.lean`).
5. **Scalar enclosures:** the twisted trace tables with their sector and filter lemmas (spin-1); the
   Ising case is done in `Ising/Seed.lean`.

Estimates: $900\text{--}1{,}950$ Lean lines for the Ising point, $8{,}000\text{--}15{,}000$ for the
spin-1 route ([`METHOD.md`](METHOD.md)). Target Design C first (real arithmetic only, largest ring
12).

## Reproduce

Requires Julia 1.12 (`export OPENBLAS_NUM_THREADS=1`) and `lake` (via `elan`) for the pinned Lean
toolchain. Run from `models/lean_seeds/`; wall times below were measured on an Apple M5 laptop:

```sh
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. test/runtests.jl         # 30 s: 1,582 unit tests across all five studies (Julia only)
julia --project=. check.jl                 # 1 min: exact rational replays + data checks, then
                                           #   `lake build SeedKernel` and the Lean axiom audit
julia --project=. check.jl --full          # 4.7 min: also regenerate the spin-1 Lean witnesses,
                                           #   rebuild SeedKernel from clean, and run all benchmarks
julia --project=. certify.jl [study ...]   # 75 s: regenerate exact JSON records and Lean witnesses
```

Running `certify.jl` executes `seed_design` (5 s), `anchor` (3 s), `ising` (9 s), and
`spin1_kernel` (55 s; run last because it records the SHA-256 hashes of shared driver files).
`check.jl` uses the Mathlib build cache in `../../lean/` (run `lake exe cache get` inside `lean/`
if not yet initialized); with Mathlib cached, `lake build SeedKernel` from clean takes about 80 s,
the axiom audit (`lake env lean SeedKernel/Audit.lean`) takes 5–30 s, and the `SeedKernelBench`
files take about 1.7 min.

Exploratory and table-generation scripts by study (writing to `results/<study>/`, with temporary
caches in the gitignored `cache/`, `upstream/`, `.cache/`, and `runtime/` directories):

```sh
julia scripts/fetch_paper.jl         # fetch upstream paper TeX -> upstream/ (network; tools/fetch_upstream.jl)
# seed_design: spectra (70 s for n <= 10), parameter searches, diagnostics, trial contraction, summary
julia --project=. scripts/seed_design/spectra.jl 4 5 6 7 8 9 10
julia --project=. scripts/seed_design/search_existing.jl   # writes results/seed_design/raw/
julia --project=. scripts/seed_design/search_new.jl 8      # or 10
julia --project=. scripts/seed_design/refine.jl
julia --project=. scripts/seed_design/diagnostics.jl
julia --project=. scripts/seed_design/shift_search.jl      # 8 s
julia --project=. scripts/seed_design/trial_scan.jl        # network
julia --project=. scripts/seed_design/summarize.jl         # requires raw/, spectra, upstream/
# spectrum: 21 rings (2.6 min), analysis, fixed-point filters, packing counts, cost models (~10 s)
julia --project=. scripts/spectrum/spectra.jl
julia --project=. scripts/spectrum/analyze.jl
julia --project=. scripts/spectrum/filters.jl float        # also fixed, sweep
julia --project=. scripts/spectrum/packing.jl
julia --project=. scripts/spectrum/anchor_probe.jl
julia --project=. scripts/spectrum/summarize.jl
julia --project=. scripts/spectrum/runtime_estimates.jl
julia --project=. scripts/spectrum/provenance.jl
# anchor: parameter scans, ring spectra, analysis, tail bins, homotopy continuation (~40 min), tables
julia --project=. scripts/anchor/scan.jl
julia --project=. scripts/anchor/count_scan.jl
julia --project=. scripts/anchor/spectra.jl 8 10 12
julia --project=. scripts/anchor/analyze.jl
julia --project=. scripts/anchor/tail_bins.jl
ANCHOR_FRESH=1 julia --project=. scripts/anchor/homotopy.jl 12 I binned   # per n, g
julia --project=. scripts/anchor/homotopy_sectors.jl      # n = 8, 10; then `I`, `P` for n = 12
julia --project=. scripts/anchor/report_tables.jl         # updates METHOD.md §3 tables and homotopy_summary.csv
# spin1_kernel: Float64 structure/anchor tables (5 s), Lean kernel benchmarks, Mathlib audit
julia --project=. scripts/spin1_kernel/experiments.jl
julia --startup-file=no scripts/spin1_kernel/benchmark.jl 300   # ~2.5 min
julia scripts/audit_mathlib.jl       # Mathlib audit; searches lean/.lake/packages/mathlib, needs rg
# ising: two-level feasibility sandwich, cost model, 43-row table, all-model inventory
julia --project=. scripts/ising/feasibility.jl
julia --project=. scripts/ising/costs.jl --reuse           # omit --reuse to re-run kernel microbenchmarks
julia --project=. scripts/ising/report_tables.jl
julia --project=. scripts/ising/inventory.jl               # see "Layout" for required spectrum files
```

## What is rigorous

- **Lean 4 kernel (`decide +kernel`, standard axioms only):** All 40 declarations audited in
  `lean/SeedKernel/Audit.lean` (verified automatically by `check.jl`) depend at most on standard
  Lean axioms (`propext`, `Classical.choice`, `Quot.sound`), with no `sorryAx` and no
  `Lean.ofReduceBool` (`native_decide`). The computational witness checks (spin-1 $H_6, H_8, H_{10}$
  moment inequalities and tests, and Ising Gram, residual, and classical anchor checks) depend only
  on `propext` or no axioms at all. Although `Ising/Seed` imports `IsingLongitudinal` (whose `seeds`
  declaration contains a `sorry`), the axiom audit confirms that no audited declaration depends on
  it.
- **Unconditional (Ising):** `seed_spatial`, `seed_thermal`, `first_point_seed`, and `ising_gap`
  have no hypotheses. **Conditional (spin-1):** the witnesses in `SeedKernel/Spin1/H*.lean`
  establish exact integer inequalities on packed moments, which become physical eigenvalue bounds
  only once the packing and spectral obligations above are discharged.
- **Exact rational arithmetic in Julia (not in Lean):**
  - `seed_design`: `Rational{BigInt}` filter exterior checks, integer-root and series enclosures,
    and a line-by-line transcription of Lean's `ChainCertificate.Valid`. The length-144 MPS trial
    in Design B is a stored exact Nemo contraction (`results/seed_design/trial_scan.json`), not
    recomputed by the checker.
  - `anchor`: exact `Rational{BigInt}` block $\mathrm{LDL}^\top$ decompositions of the $81 \times 81$
    four-site remainders, half-chain eigenvalue enclosures from `Float64` eigenvector proposals via
    bounded `Int128` polar decomposition and Weyl's inequality, exact tensor-sum comparison
    counts, and degree-256 rational Taylor upper bounds on the binned tail.
- **Untrusted generators:** `Float64` LAPACK proposes the Ising eigenvectors and the anchor
  eigenvectors; wrong proposals fail the exact checks. The Julia replays of kernel-checked facts
  (`check_spin1_kernel`, `check_ising`) are diagnostics.
- **Floating-point only:** all spectra, counts, cutoffs, filter residuals, searches, Pareto
  fronts, homotopy continuations, inventories, and cost/time forecasts.
- **Platform:** LAPACK proposals differ across platforms; the stored Ising and anchor records were
  regenerated on the M5 machine and reproduce bit for bit there. Other platforms may give
  different but equally valid proposals.

## Layout

| Path | Content |
|---|---|
| `certify.jl`, `check.jl` | Regenerates exact records and Lean witnesses by study; verifies them (fast, or `--full` for a complete rebuild and benchmarks) |
| `src/LeanSeeds.jl` | Top-level module `LeanSeeds`, loading one submodule per study |
| `src/Common.jl` | Shared paths, SHA-256 source hashing, $S^z$ basis utilities, Lake runner, and axiom-audit parser |
| `src/seed_design.jl` | `SeedDesign`: `Float64` search, exact rational filters/seeds/chains, `certify`, and `check` |
| `src/spectrum.jl` | `Spectrum`: twisted $`S^z`$ Hamiltonians, low-energy spectra, fixed-point Chebyshev filters, and data checks |
| `src/anchor.jl` | `AnchorStudy`: integer bonds, boundary comparison matrices, exact $\mathrm{LDL}^\top$ and half-chain enclosures, `certify` and `check` |
| `src/spin1_kernel.jl` | `Spin1Kernel`: witness generator, Lean source emitter, witness inventory, and checker |
| `src/ising.jl` | `IsingSeed`: eigenvector proposal generator, `Generated*.lean` emitter, exact Julia replay, and checker |
| `scripts/<study>/` | Exploratory, benchmark, and table-generation drivers listed under "Reproduce" |
| `scripts/fetch_paper.jl` | Shared: runs `tools/fetch_upstream.jl --paper` (paper sources into `upstream/`) |
| `scripts/audit_mathlib.jl` | `spin1_kernel` Mathlib audit; kept outside `scripts/spin1_kernel/`, whose `.jl` files are hash-pinned by `results/spin1_kernel/witnesses.json` |
| `test/runtests.jl` | Master test suite running `test/<study>.jl` (one module each) |
| `results/seed_design/` | `certificate.json`, `compact_witnesses.json` (exact); `trial_scan.json`, `pareto.csv`, `summary.json`, `diagnostics.json`, `shift_search.json` (search records) |
| `results/spectrum/` | `spectrum_<n>_<twist>.json`, `analysis.json`, `summary.json`, CSV tables, kernel logs, `provenance.json` |
| `results/anchor/` | `certificate.json` (exact); `optimized.json`, `count_optimized.json`, `binned_tails.json`, `sectors.csv`, `homotopy_summary.csv` |
| `results/spin1_kernel/` | `witnesses.json`, `exact_witnesses.tsv`, `axioms.txt`, `benchmarks*.tsv`, `mathlib-audit.txt`, `paper-source-sha256.txt` |
| `results/ising/` | `certificate.json` (integer proposals); `feasibility.json`, `inventory.{json,tsv}`, `costs.json` |
| `../../lean/SeedKernel/` | Lean 4 `SeedKernel` and `SeedKernelBench` libraries (file map in [`lean/README.md`](../../lean/README.md)) |

**Hash-pinning notes:** `check.jl` verifies the source hashes recorded in
`results/spin1_kernel/witnesses.json` against the current `lean/SeedKernel/Spin1/` files, Lake
configuration, and shared Julia drivers (`certify.jl`, `check.jl`, `src/Common.jl`,
`src/spin1_kernel.jl`), and checks `results/ising/certificate.json` against `src/ising.jl`. If you
modify any of those code files, rerun `certify.jl spin1_kernel` (and `certify.jl ising`). No
Markdown file (`.md`) is hash-pinned. The script `scripts/ising/inventory.jl` reads stored spectrum
files across other `models/` directories; two large spectrum archives
(`heisenberg_ladder/results/verified_L8.json.gz` and `spin1_dimerized/results/spectrum_d9_c4.json`)
are gitignored and must be regenerated in their respective model directories before re-running
`inventory.jl` from scratch.

## History

- The five studies were originally developed in separate directories, Julia projects, and PRs:
  PR #24 (`seed_design`), PR #22 (`spectrum`), PR #27 (`anchor`), PR #21 (`spin1_kernel`, then
  `lean_prototype`), and PR #25 (`ising`, then `models/lean_seed_small`), with two standalone Lake
  projects for the Lean prototypes.
- Consolidation: one Julia project here (the union of the environments, same package versions);
  one module per study with shared paths, $`S^z`$ basis, Lake driver, and source hashing; the paper
  fetchers merged into `scripts/fetch_paper.jl`. The two Lake projects became the opt-in library
  `SeedKernel` of `lean/`, with one shared packed-lane module (`SeedKernel/Packed.lean`, now also
  under the Ising `Packed` checks) and `SeedBridge` importing `CertifiedGaps` directly. All exact
  records were regenerated: every exact number is identical (only provenance hashes, timings, and
  the `seed_design` trial-file paths changed), and the Lean witness files are byte-identical apart
  from module names.
- Layout tidy: `scripts/fetch_paper.jl` now forwards to the shared
  `tools/fetch_upstream.jl --paper` (one place for the pinned commit; it also skips files that are
  already present with the pinned hash), and `scripts/spin1_kernel/audit_mathlib.sh` was ported to
  `scripts/audit_mathlib.jl` (byte-identical output). The port sits outside `scripts/spin1_kernel/`
  so that the `spin1_kernel` source-hash set is unchanged. `scripts/spectrum/KernelBench.lean` stays
  where it is: `scripts/spectrum/provenance.jl` hashes it into `results/spectrum/provenance.json`,
  and the Lean benchmarks under `lean/` are not edited from model directories.
- The Ising seed closed: `Ising/Spectral`, `Ising/Hamiltonian`, `Ising/Bridge`, and `Ising/Seed`
  prove `first_point_seed` and the unconditional `ising_gap` from the existing kernel checks (no new
  witnesses); `SeedBridge`, `Arithmetic`, and `ResidualLemma` were retired.
