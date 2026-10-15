# Toward a Lean-checked seed: route memos

This document collects the five route memos written prior to computation for the studies in
`models/lean_seeds/` (`seed_design`, `spectrum`, `anchor`, `spin1_kernel`, and `ising`), followed
by a summary of their outcomes. Across all five studies, the guiding architecture is the same:
Julia acts as an untrusted numerical search tool and witness generator, while the ultimate target
is a self-contained proof inside the Lean 4 kernel (`decide +kernel`, with no `native_decide` and
no custom axioms). Floating-point calculations below are exploratory and non-certifying; exact
rational and integer checks establish finite algebraic identities that become unconditional gap
theorems only when paired with their corresponding Lean bridge lemmas.

## 1. `seed_design`: which finite inputs does the spin-1 seed need?

1. **Even-length real-twist restriction (`J/N` split):** Restrict existing source-1 thermal data
   to even lengths up to $n = 10$ or $n = 8$. Optimize nonnegative polynomial moment filters
   $F(x) = x^{2s} q(x)^2$ (with even $q$) and preserve the threefold $N$-sector multiplicity using
   only the untwisted ($I$) and $\pi$-twisted ($P$) traces via $J = (Z + 3P)/4$ and
   $N = (Z - P)/4$. This avoids complex $2\pi/3$-twist ($C$) arithmetic altogether. Because kernel
   verification cost scales as $\sum 3^n K_n$ (where $K_n$ is the low-energy state count), smaller
   maximum ring size $n_{\max}$ and lower polynomial degree reduce both the number of trace tables
   and scalar proof complexity. *Risk:* the modulus cap on the combined sector $J = I + 2O$ may
   remain too far above $1$, causing its high power at the MPS trial length $N$ to blow up the
   defect bound.
2. **Full three-twist decomposition (`I/O/N` split):** Compare the full three-twist decomposition
   (using $I, P,$ and $C$ to separate $I, O,$ and $N$) against the two-twist $J/N$ split and a
   twist-free bound. Isolating the $O$ sector at shorter chain lengths ($n \le 10$) may tighten the
   $I$-sector modulus cap enough to more than pay for the complex Eisenstein arithmetic of the $C$
   twist.
3. **Temperature lifting:** Raise the inverse temperature $\beta_0 = r \beta_{\mathrm{table}}$
   before extracting seed purities. Scan over base inverse temperatures $\beta_{\mathrm{table}}$,
   lifting multipliers $r$, intermediate moment orders $k$, and MPS trial lengths $N$. A longer
   trial ring improves spatial isolation of the ground state, at the cost of raising an imperfect
   upper modulus cap to a higher power in the thermal purity estimate.
4. **Per-site energy shift $\kappa$ ($a$):** Vary the energy shift within the margin permitted by
   the variational MPS upper bound. Because shifts rescale every transfer eigenvalue and moment
   uniformly, their only effect comes from trading variational margin against trace-table error
   scaling.
5. **Relaxed table tolerances and uncapped sectors:** Compare capped versus uncapped sector mass
   estimates, lower-degree polynomial filters, and relaxed absolute trace tolerances $\epsilon$.
   Accepting a slightly smaller certified gap lower bound is well worthwhile if it eliminates the
   $n = 12$ trace tables ($3^{12} = 531{,}441$ states) or widens the allowable trace error by an
   order of magnitude.

**Chosen workflow:** Reuse the stored exact thermal trace centers first, compute new small-ring
spectra by `Float64` magnetization-block exact diagonalization or Lanczos only where needed to
explore other temperatures, search for simplified scalar designs, and replay the best 2–3 designs
in exact rational arithmetic down to `ChainCertificate.Valid`.

## 2. `spectrum`: how many low states, how precise, how expensive?

1. **Dense diagonalization through $n = 8$:** Full dense diagonalization in each $S^z$ sector
   through $n = 8$ provides an independent multiplicity check and complete thermal sums at cost
   $O(d_m^3)$ per sector.
2. **Sparse Lanczos/Arnoldi through $n = 12$:** Restarted Krylov eigensolvers through $n = 12$
   compute the low-energy spectrum above the tail cutoff at cost proportional to $\mathrm{nnz}(H)$
   times Krylov iterations, validated against dense spectra at $n \le 8$ and checked with guard
   eigenvalues above every counting threshold.
3. **Fixed-point Chebyshev block filtering:** Test the proposed Chebyshev polynomial filtering
   recipe starting from deterministic computational-basis product states, incorporating subspace
   oversampling and periodic Rayleigh–Ritz orthogonalization to prevent rank collapse. Measure
   total matrix-vector ($Hv$) products across cycles, and compare 48-bit (Q48) and 112-bit (Q112)
   floor-rounded integer recurrences against `Float64`.
4. **Tail cutoff and tolerance allocation:** Quantify how the absolute shifted-trace tolerance
   $\delta$ translates into required individual eigenvalue accuracy $\varepsilon$ across chain
   lengths $n = 4\text{--}12$, and compare sectorwise tail cutoffs against the uniform
   dimension-times-exponential bound $(3^n - k)e^{-\beta c}$.
5. **Full base-3 vs. compressed $S^z$ bit-packing:** Full ternary indexing ($3^n$ lanes) implements
   nearest-neighbor spin flips as fixed bit shifts across $2n$ masks, whereas compressed $S^z$
   indexing reduces lane count (from $531{,}441$ to $73{,}789$ at $n = 12, S^z = 0$) at the price
   of irregular rank permutations. Count the exact number of displacement masks required in the
   compressed basis rather than assuming both small lane count and $O(n)$ shifts simultaneously.

**Chosen route:** Compute dense and sparse reference spectra across all 21 ring/twist pairs, run
Q48 and Q112 Chebyshev filtering experiments across all active sectors, count exact displacement
masks in both basis layouts, and benchmark bulk `Nat` operations in the Lean kernel to build a
calibrated runtime model.

## 3. `anchor`: can a cheap comparison anchor the eigenvalue count?

1. **Bare open half-chains:** Cutting two opposite bonds of the ring and replacing each cut bond by
   its minimum eigenvalue $-2I$ yields $H \ge H_{\mathrm{open}, l} \otimes I + I \otimes H_{\mathrm{open}, l} - 4I$
   ($l = n/2$). Because the right-hand side is a Kronecker sum, its spectrum is the pairwise sum
   of half-chain eigenvalues ($S^z$ blocks of dimension at most $141$ at $l = 6$), at the cost of a
   $\sim 2$-unit surface energy penalty that inflates the state count below the thermal cutoff.
2. **Tuned two-site $\mathrm{SU}(2)$-invariant boundary corrections:** Parameterize boundary
   operators $B = u h + v h^2$ at the ends of each half-chain so that the four-site cut remainder
   $R = h_{23} - B_{12} - B_{34} - cI$ is positive semidefinite on $\mathbb{C}^{81}$. Because $R$
   decomposes into $S^z$ blocks of size at most $19 \times 19$, $R \ge 0$ can be certified by an
   exact rational $\mathrm{LDL}^\top$ factorization. Optimize $(u, v, c)$ directly to minimize the
   number of comparison eigenvalues below the thermal tail cutoff $\rho$, rather than maximizing
   the ground-state lower bound.
3. **Three or more disjoint open blocks:** Cutting the ring into three or four pieces reduces the
   subgraph dimensions further, but introduces additional cut boundaries that may worsen the total
   surface energy loss and state overcount. Measure this tradeoff directly.
4. **Overlapping Anderson clusters:** Translating a short open cluster around the ring yields
   tighter scalar lower bounds on the ground-state energy $E_0$, but because overlapping cluster
   Hamiltonians do not commute, their sum does not have a simple tensor-sum excited spectrum to
   supply a counting anchor $\lambda_{k+1}(H) \ge \rho$.
5. **Effective edge-spin projections:** Projecting onto low-energy edge states yields variational
   upper bounds, whereas a lower comparison operator requires certified control of the discarded
   subspace and cross-couplings.
6. **Homotopy continuation from $H_0$ to $H$:** When the comparison operator $H_0 \le H$ has
   $K > k$ eigenvalues below the target cutoff $\rho$, a one-shot comparison count only bounds
   $\lambda_{K+1}(H) \ge \rho$. Interpolating along $H(t) = H_0 + t(H - H_0)$ for $t \in [0, 1]$
   uses operator monotonicity to propagate upper counting barriers while inverse-Ritz lower bounds
   shed levels that cross above $\rho$.

**Chosen route:** Implement Routes 1 and 2 with exact rational $\mathrm{LDL}^\top$ certificates and
verified half-chain spectral enclosures, tighten the tail bound by binning the comparison density
of states, compare against multi-block partitions (Route 3), and test adaptive homotopy
continuation (Route 6) across all active sectors at $n = 8, 10, 12$.

## 4. `spin1_kernel`: does packed arithmetic fit the kernel?

1. **Full base-3 radix packing (`Packed.lean`):** Pack all $3^n$ amplitudes into a single
   arbitrary-precision `Nat` with radix $B = 2^w$. Digit-pattern masks have closed-form geometric
   series representations, nearest-neighbor exchange hops become fixed bit shifts, and inner
   products reduce to a single large integer multiplication via butterfly lane reversal.
2. **Compact magnetization-sector packing:** Audit the exact number of rank-difference shift
   classes in the $S^z = 0$ sector for $n = 6, 8, 10, 12$ to test whether compact packing is
   viable without a custom routing network.
3. **Marshall-transformed positive power filter:** On even untwisted rings, the Marshall sign
   transformation makes all off-diagonal entries of $nI - H_{\mathrm{M}}$ nonnegative, allowing
   unsigned power-filter trial generation $v = (nI - H_{\mathrm{M}})^q e_{\text{N\'eel}}$. For
   general signed and twisted applications, represent signed vectors by separate positive and
   negative `Nat` limbs with explicit carry and floor-division bounds.
4. **Exact integer Rayleigh–Ritz and residual witnesses:** Kernel-check the exact integer moment
   inequalities $G > 0$, $\mathrm{den}\cdot A \ge (\mathrm{den}\cdot n + \mathrm{num})G$, and
   $\mathrm{epsDen}^2 (CG - A^2) \le \mathrm{epsNum}^2 G^2$ for $H_6, H_8, H_{10}$, and test
   whether the untuned cut anchor satisfies the strict separation condition $\theta_k + \varepsilon < \rho$
   required by the Weyl counting lemma.
5. **Mathlib spectral API audit and formalization budget:** Audit the pinned Mathlib commit for
   Hermitian spectral theorem, Courant–Fischer min–max, Weyl perturbation, Cauchy interlacing,
   and Jensen/Peierls inequalities to establish an accurate line-of-code and runtime budget.

**Chosen route:** Build the executable `Std`-only Lean prototype for full base-3 packing, benchmark
$n = 6, 8, 10, 12$ matrix-vector and dot products in `decide +kernel` (including hoisted masks and
division-free reversal), kernel-check $H_6, H_8, H_{10}$ moment witnesses, and audit Mathlib.

## 5. `ising`: the cheapest model end to end

1. **Full dense eigen-certificate:** For the longitudinal Ising chain ($n_0 = 6, 2n_0 = 12$), the
   full Hilbert spaces have dimensions $64$ and $4096$. Certifying all $4096$ eigenvectors in the
   Lean kernel avoids needing a tail cutoff but requires quadratic storage and cubic Gram checking.
2. **Low-energy Ritz vectors plus a classical diagonal anchor:** Bound the off-diagonal transverse
   field operator-wise by $-g \sum_j X_j \ge -n g I$, reducing the counting anchor to an exact
   combinatorial count of classical Ising spin configurations below the domain-wall threshold. At
   weak transverse field ($g = 1/4, h = 1/20, \beta_0 = 97/4$), only the two ferromagnetic ground
   states lie below the first domain-wall excitation energy $\rho_n = (82 - 26n)/20$, so only two
   trial vectors per ring are required.
3. **Characteristic polynomial / Sturm chain on the small ring:** Fraction-free characteristic
   polynomials could certify $n = 6$ ($d = 64$), but remain prohibitive at $2n = 12$ ($d = 4096$).
4. **Perturbation from the classical Ising limit:** Pure operator perturbation bounds without
   variational eigenvectors lose the delicate finite-size cancellation between $2\log Z_6(\beta_0)$
   and $\log Z_{12}(\beta_0)$ required for the spatial purity seed.
5. **Partition-function enclosure interface:** Formalize the rational arithmetic step from explicit
   $\log Z$ interval hypotheses to the exact `IsingLongitudinal.seeds` conjunction in Lean, cleanly
   separating the kernel-checked linear algebra witnesses from the remaining analytical spectral and
   exponential lemmas. This is partial progress only: the interface hypotheses must remain visible
   until their matrix-exponential semantics are proved.

**Chosen route:** Combine Route 2 (two-vector integer residual + Gram + classical domain-wall anchor
checked by `decide +kernel` at $n = 6$ and $n = 12$, in both packed and chunked forms) with Route 5
(the partition-function enclosure interface), and survey the low-energy state counts $k$ across all
43 certified parameter points in the repository to quantify where the low-rank Ritz + anchor
strategy scales.

## Outcome

### `seed_design`

Route 2 (the three-twist $I/O/N$ split using the $C$ twist at $n = 6, 8, 10$) and Route 1 at
$n \le 12$ (even lengths $4, 6, 8, 10, 12$ with real twists $I, P$ only) both succeeded and were
verified in exact rational arithmetic (`results/seed_design/certificate.json`). Design A uses just
9 traces at $n \in \{6, 8, 10\}$ with absolute tolerance $3 \times 10^{-4}$ to certify
$\Delta_L > 0.005966152$ for even $L \ge 60$ (an $8\times$ reduction in total Hilbert-space
dimension relative to source 1), while Design C uses 10 real-twist traces at even $n \le 12$ with
tolerance $10^{-3}$ to certify $\Delta_L > 0.039036761$ for even $L \ge 60$. Conversely, dropping
$n = 10$ ($n_{\max} = 8$) or dropping the $C$ twist at $n_{\max} = 10$ fails across the entire
tested filter and temperature-lifting grid, and shifting $\kappa$ ($a$) merely rescales the error
budget without reducing work.

### `spectrum`

Across all 21 ring/twist pairs, the number of physical states below the thermal tail cutoff is
small ($26$ for $I$ and $34$ for $P$ at $n = 12, \beta = 49/4$; $22, 30, 26$ for $I, P, C$ at
$n = 10, \beta = 21/2$), and degree-16 Q48 integer Chebyshev filtering meets the required
eigenvalue residual tolerances in every active sector. However, the study exposed two critical
bottlenecks in the naive kernel plan: (1) compressed $S^z$ indexing requires thousands of
displacement masks ($3{,}816$ at $n = 12, S^z = 0$; $13{,}620$ at $S^z = 2$) rather than $\sim 50$,
making full $3^n$ ternary packing with $2n$ shifts preferable; and (2) the bare two-half-chain
comparison leaves $242$ states in $S^z = 0$ below the $n = 12$ tail cutoff compared with $10$
physical states.

### `anchor`

Optimizing the two-site $\mathrm{SU}(2)$-invariant boundary terms $B = u h + v h^2$ directly for
the cutoff state count (Route 2) and verifying the $81 \times 81$ four-site remainder by exact
rational $\mathrm{LDL}^\top$ decomposition cuts the $n = 12, S^z = 0$ comparison count from $191$
(bare cuts) and $212$ (ground-energy optimized) down to $121$. Replacing the uniform dimension tail
by a 32-bin histogram of the comparison spectrum reduces the initial $S^z = 0$ count further to
$49$ ($183$ across all $S^z$ sectors). Because $49$ still exceeds the $8$ physical $S^z = 0$
levels, a one-shot comparison count is obstructed by the missing-level counterexample, making
homotopy continuation (Route 6) or a complementary subspace certificate necessary. Numerical
homotopy continuation succeeded in all 26 active sectors ($n = 8, 10, 12$), whereas cutting into
three or more blocks (Route 3) worsened the counts.

### `spin1_kernel`

Full base-3 radix packing (Route 1) works cleanly inside the Lean kernel: dense $n = 12$
matrix-vector and dot products check via `decide +kernel` in $9.6\text{--}26\text{ s}$ on an
Apple M5, and exact power-filter moment witnesses for $H_6, H_8, H_{10}$ compile in seconds.
Hoisting digit masks and replacing big-integer divisions by shift-and-OR bit reversal (`Optimized.lean`)
sped up $n = 12$ dot products by $5.7\times$. However, testing the untuned cut anchor against the
exact $H_6$ and $H_8$ spectra confirmed that $\max_k (\lambda_{k+1}(H_0) - \lambda_k(H)) < -1.38$,
so the separation premise of the one-shot Weyl counting lemma fails without a tuned anchor or
continuation.

### `ising`

Combining Route 2 (two approximate eigenvectors plus the classical domain-wall anchor, with only two
levels per ring) and Route 5 (the partition-enclosure interface) at
$(g, h) = (1/4, 1/20), (n_0, \beta_0) = (6, 97/4)$ produced a complete, kernel-checked arithmetic
verification of the integer Gram, residual, and anchor inequalities at both $n = 6$ and $n = 12$ (in
both packed-lane and 64-chunk direct forms). Later the seed was closed without new witnesses: the
matrix identification, counting by comparison (instead of min–max), and elementary `exp` bounds give
the unconditional `ising_gap` (README, "ising"; `lean/SeedKernel/Ising/Seed.lean`). The 43-point
cross-model inventory showed that while this low-rank approach is fast for the 30 longitudinal Ising
points and the $n_0 = 4$ spin-1 dimerized point, full $k \times k$ Gram checks reach hours to
$\sim 1{,}000$ hours for models with $65{,}536$- or $531{,}441$-dimensional seeds at smaller
$\beta_0$ ($300\text{--}2{,}465$ states below the tail cutoff), so Routes 1–2 do not scale there
without symmetry block-diagonalization or cluster methods.

### Overall synthesis

The generator-and-kernel architecture is viable: Julia generates witness vectors and rational
certificates, and packed `Nat` arithmetic or chunked evaluation checks them inside Lean's kernel
with no new axioms. One seed is closed: the longitudinal Ising point $(g, h) = (1/4, 1/20)$ now has
an unconditional, kernel-checked gap theorem. For the spin-1 Heisenberg chain no seed is closed;
Design C of `seed_design` (even $n \le 12$, real twists $I, P$) is the recommended target. The
bottleneck there is the spectral counting anchor, for which the local PSD comparison of `anchor` is
exact and cheap but leaves too many levels for a one-shot count, so a certified continuation (or a
better comparison) is the next research question. The missing Lean linear algebra (min–max, Weyl,
interlacing) and the packing/Hamiltonian bridges are the formalization cost (README, "Open
obligations").
