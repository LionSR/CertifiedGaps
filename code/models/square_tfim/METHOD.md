# Uniform gaps for the square-lattice transverse-field Ising model

## Hamiltonian and small periods

For even periods $n,m\ge2$, with both indices read periodically, the Pauli Hamiltonian is

```math
H_{n,m}=-\sum_{x=0}^{n-1}\sum_{y=0}^{m-1}
 (Z_{x,y}Z_{x+1,y}+Z_{x,y}Z_{x,y+1}+gX_{x,y}).
```

Thus an undirected bond occurs twice in a direction of period two. This convention is
essential: it is the trace of the same transfer matrix at periods two and four. For periods
at least four the expression is the usual square-lattice nearest-neighbor Hamiltonian.
The conclusion concerns square tori with both periods tending to infinity.

The simultaneous Hadamard rotation exchanges X and Z. In this frame, each row has
dimension $2^m$, its on-site term is $-\sum_y X_yX_{y+1}-g\sum_y Z_y$, and its inter-row
bond is $-\sum_y X_y\otimes X_y$. The factors $C_y=X_y$ are exact integer matrices.
`SquareIsing.row_rp` verifies this identity using the shared core. Transposition of the
square lattice gives the same construction in the other direction.

## Spatial inequalities without a Trotter error assumption

Put $\delta=\beta/M$. In the Z basis, the product formula

```math
Z_M(n,m,\beta)=\mathrm{Tr}
 \left(e^{\delta\sum_{\langle u,v\rangle}Z_uZ_v}
       e^{\delta g\sum_u X_u}\right)^M
```

is a classical partition sum over periodic spin histories. The spatial edge weights are
$e^{\delta st}$, while each temporal edge has weight
$w(s,t)=\langle s|e^{\delta gX}|t\rangle$, equal to $\cosh(\delta g)$ or
$\sinh(\delta g)$. All weights are positive when $g>0$.

For a complete column history $s$, let $D(s)>0$ be the product of weights internal to that
column. The interaction between adjacent histories is

```math
K(s,t)=\exp\left(\delta\sum_{y,k}s_{y,k}t_{y,k}\right),\qquad
T_x=D^{1/2}KD^{1/2},\qquad Z_M(n,m,\beta)=\mathrm{Tr}\,T_x^n.
```

The matrix K is a tensor product of matrices with diagonal $e^\delta$ and off-diagonal
$e^{-\delta}$, hence is positive semidefinite. The matrix $T_x$ is real symmetric and
independent of n. The analogous statement holds in y. In particular the normalized even
powers give spatial purities in $(0,1]$, their squaring inequality, and the monotonicity of
spectral power means used for interpolation between even periods.

At each fixed finite torus and inverse temperature, the finite-dimensional Lie product
formula gives $Z_M\to Z=\mathrm{Tr}\,e^{-\beta H}$. Every spatial inequality just listed
involves finitely many strictly positive partition functions, so it passes to this limit.
This proves the scalar spatial assumptions needed in the cube argument. A limiting transfer
Hilbert space need not be constructed. Temporal squaring is applied directly to the exact
Gibbs operator, after taking the limit; no temporal squaring identity for a fixed M is assumed.

## The three seed defects

The finite spectral calculation encloses four numbers:

```math
a=\log Z(2,2,\beta),\quad b=\log Z(4,2,\beta),\quad
c=\log Z(4,4,\beta),\quad d=\log Z(4,4,2\beta).
```

The upper bounds on the defects are obtained from $1-e^{b-2a}$, $1-e^{c-2b}$, and
$1-e^{d-2c}$, with outward interval rounding. These quantities are invariant under adding
a constant energy per site.

The core `VerifiedSpectra.ring_blocks` splits each row-ring Hamiltonian by global spin parity,
longitudinal momentum, and commuting transverse reflections. For width four these are the
exchanges of opposite vertices 0 with 2 and 1 with 3. The lattice reflection makes the
momentum blocks real. Every block is enclosed by `verified_eigen`; all multiplicities are
included and must sum to $2^{nm}$. The log partition sums use the core interval routine.
This is a new two-dimensional certificate, not an invocation of the one-dimensional gap
theorem: only the core model, reflection-positivity and spectral arithmetic are reused.

## Exact recurrence and all even square periods

Write $e=(e_x,e_y,e_\tau)$ and $f(z)=z^2/[2(1-z)^2]$. The update is

```math
a=f(1-(1-e_x)^4(1-e_y)^2(1-e_\tau)),
\qquad b=f(1-(1-e_y)^4(1-e_\tau)^2(1-a)),
\qquad c=f(1-(1-e_\tau)^4(1-a)^2(1-b)).
```

These are exactly the staggered updates in the repository's 2+1 dimensional note. Upward
rational rounding after each component preserves the bounds, since all expressions are
increasing on the admissible domain. The three new defects occupy the same staggered
positions after doubling both spatial periods and inverse temperature.

At scale j, set $\beta_j=2^j\beta$ and
$Q_j=(1-e_{x,j})^4(1-e_{y,j})^2(1-e_{\tau,j})$.
For each even square side $2^{j+1}\le L\le2^{j+2}$, spatial moment interpolation gives
$P_\tau(L,L,\beta_j)\ge Q_j$. If $Q_j>1/2$, the largest Gibbs weight is at least $Q_j$,
the ground state is unique, and

```math
e^{-\beta_j\Delta_{L,L}}\le\frac{1-Q_j}{Q_j},\qquad
\Delta_{L,L}\ge\frac{\log(Q_j/(1-Q_j))}{\beta_j}.
```

Once $r=\max e_j\le1/1000$, the note gives $r'\le25r^2$. Indeed the coefficient bound
reduces to the exact integer inequality $49{,}000{,}000<50\cdot993^2$. Starting at that
scale, $r_k\le(25r)^{2^k}/25$, and intermediate-period interpolation gives defect at most
$7r_k$. Since $7/[25(1-7r_k)]<1$, every subsequent window has

```math
\Delta_{L,L}>\frac{\log(1/(25r))}{\beta_j}.
```

The minimum of the finitely many earlier window bounds and this tail bound covers every
even square side at least two. All logarithms are bounded below by the core exact rational
atanh-series routine `ExactArithmetic.log_enclosure`. In particular, using initial defects
larger than $10^{-3}$ does not assume the small-defect theorem outside its stated domain.

## Refinement and finite-volume energies

Replacing the three defects by their maximum loses information. The refined calculation
retains the three-component recurrence for twelve steps and then applies the same scalar
tail. Upward rounding uses denominator $10^{2000}$; the original denominator $10^{80}$
would dominate the defects before twelve steps. A numerical search over
$\beta g=3,3+1/40,\ldots,6$ selects the temperatures below. Fresh interval partition sums
certify these choices independently of the search.

| $g/J$ | $\beta_0J$ | original bound $/J$ | refined bound $/J$, rounded down |
|---:|---:|---:|---:|
| 18 | $53/240$ | 0.5943744623 | 0.7437803188 |
| 20 | $157/800$ | 2.1736330392 | 2.5137666679 |
| 24 | $79/480$ | 4.8055785105 | 5.5217326453 |

The original seed, recurrence, and spectral records are retained exactly. Certificate
format v2 adds a `refinement` record and rational enclosures of the ground energy, energy
per site, and gap for each finite torus. Both chains and every new enclosure are recomputed
by the checker. The finite-volume values in the README come from `ground_energy_enclosure`
and `spectral_gap_enclosure` applied to the complete verified spectra.

The bootstrap bounds remain well below the sixteen-spin gaps. At $g/J=24$, for example,
the finite gap exceeds $43.91119J$, while the uniform bound is $5.52173J$. Spectral enclosure
errors are below $10^{-6}J$ per eigenvalue. Sharper eigenvalue arithmetic alone would thus
have little effect. Larger spatial seeds or sharper directional estimates are needed for
a substantial further increase.

## General theorem and Lean formalization

Let $Z_{x,y,t}>0$ be partition data indexed by nonnegative integers counting doublings in
two spatial directions and inverse temperature. Define

```math
p_x(x,y,t)=\frac{Z_{x+1,y,t}}{Z_{x,y,t}^2},\qquad
p_y(x,y,t)=\frac{Z_{x,y+1,t}}{Z_{x,y,t}^2},\qquad
p_t(x,y,t)=\frac{Z_{x,y,t+1}}{Z_{x,y,t}^2}.
```

Assume that each ratio is at most one and satisfies the single-direction squaring
inequality: $1-p_i(v+e_i)\le f(1-p_i(v))$. These are the fields of the Lean structure
`CertifiedGaps.Spacetime.CubeData`. Cube coupling follows from cancellation of the partition
functions; it is proved from these definitions.

For a family of finite spectra $E_L$, a positive seed period $n_0$, and $\beta_0>0$, assume
the spatial interpolation inequality

```math
P_\tau(E_L,2^j\beta_0)\ge Q_j
 :=p_x(j,j,j)^4p_y(j+1,j,j)^2p_t(j+1,j+1,j)
```

for every even $L$ with $n_0 2^j\le L\le2n_0 2^j$. This is `WindowInput`; it is an
explicit hypothesis relating the partition data to the physical spectra.

**Theorem.** If $0<\varepsilon\le10^{-3}$ and

```math
p_x(0,0,0),\ p_y(1,0,0),\ p_t(1,1,0)\ge1-\varepsilon,
```

then every even $L\ge n_0$ has a unique lowest energy and

```math
\Delta_L>\frac{\log(1/(25\varepsilon))}{\beta_0}.
```

**Proof.** Cube coupling and directional squaring give the common next-scale bound
$25\varepsilon^2$. Induction yields
$u_j=(25\varepsilon)^{2^j}/25$, and Bernoulli's inequality gives $Q_j\ge1-7u_j$.
Choose the window containing $L$. For $c=\log(1/(25\varepsilon))/\beta_0$ one has
$e^{2^j\beta_0c}7u_j=7/25<1-7u_j$. The purity-to-gap inequality proves the claim.

This proof is `Spacetime.gapGT_of_cube_seeds` in
[`CubeCriterion.lean`](../../lean/CertifiedGaps/CubeCriterion.lean).
Its only axioms are Lean's standard `propext`, `Classical.choice`, and `Quot.sound`.
The formal theorem is independent of Ising and of a local Hilbert-space dimension.
It treats two spatial directions. The model's Trotter-limit construction, the derivation
of `WindowInput` from spatial transfer moments, the asymmetric finite prefix, and the
numerical seed arithmetic have not been formalized in this file.

## Moderate-field numerical checks and QMC context

The script `moderate_fields.jl` diagonalizes the complete 2 by 2, 4 by 2, and 4 by 4
spectra in floating point at $g/J=4,5,6,8$. The sixteen-spin gaps are approximately
$3.24828451J,5.44369104J,7.56254563J,11.69479751J$. These are finite-volume numerical
values, with no verified eigenvalue enclosure attached to this scan.

Across $\beta g=0.5:0.025:16$, the largest initial cube lower bounds Q are approximately
0.04075, 0.11899, 0.21982, and 0.41048. No three-direction recurrence closes at any tested
temperature. This is a negative result for the tested small seed and temperature grid.
The ground-energy doubling mismatches and full scan tables are stored in
`results/moderate_g4.json`, `moderate_g5.json`, `moderate_g6.json`, and `moderate_g8.json`.

QMC measurements of Ising excitation gaps are available in the literature.
[Schuler et al., Physical Review Letters 117, 210401 (2016)](https://arxiv.org/abs/1603.03042)
use continuous-time world-line QMC with cluster updates and extract finite-size gaps
from imaginary-time spin correlations, with system sizes up to $30\times30$ at criticality
(see their supplemental section on finite-size gap estimation). Those are literature
results; none of the numerical values computed in this directory is a QMC measurement.
The calculation here should also be distinguished from the even-parity excitation gap
shown in [Quantum phase transition dynamics (2022), Fig. 2](https://pmc.ncbi.nlm.nih.gov/articles/PMC9481121/),
which uses neural quantum states and restricts the symmetry sector.

## A graph-level high-field theorem

The small spatial seeds are unnecessary for a second, model-specific argument. The result
below follows from the ferromagnetic Ising path expansion and the FKG inequality. It is
independent of the spacetime bootstrap and does not assert a new correlation inequality.

**Theorem.** Let G be a finite nonempty graph, with symmetric couplings $J_{ij}\ge0$,
$J_{ii}=0$, and fields $g_i\ge h>0$. Set $s=\max_i\sum_jJ_{ij}$. If $h>s$, then

```math
H=-\sum_{\{i,j\}}J_{ij}Z_iZ_j-\sum_i g_iX_i
```

has a unique ground state and

```math
\Delta(H)\ge 2\sqrt{h(h-s)}.
```

The bound is uniform over graph size. Consequently the square-lattice Hamiltonian, on
any finite rectangle with open or periodic boundaries, satisfies

```math
\Delta\ge2\sqrt{g(g-4J)}\qquad(g>4J).
```

Periodic sides of length two may retain the repeated-bond convention: the aggregate
incident coupling remains $4J$. Sides of length one contribute only constant self-bonds,
which can be omitted. No even-period restriction is needed in this theorem.

### 1. A discrete imaginary-time representation

Write $V=\sum_{\{i,j\}}J_{ij}Z_iZ_j$ and $G_x=\sum_i g_iX_i$. For $a>0$ define

```math
T_a=e^{aV/2}e^{aG_x}e^{aV/2},\qquad H_a=-a^{-1}\log T_a.
```

The matrix $T_a$ is real symmetric, positive definite, and entrywise strictly positive.
Its largest eigenvalue is simple with a strictly positive eigenvector. On each fixed
finite graph, the matrix expansion at $a=0$ gives $H_a\to H$ in operator norm.

For an integer time circumference M, $\mathrm{Tr}\,T_a^M$ is a classical zero-field
ferromagnetic Ising partition function on $G\times\mathbb Z/M\mathbb Z$. The spatial
coupling is $aJ_{ij}$ and the temporal coupling at site i is
$K_i=\frac12\log\coth(ag_i)$. In particular,

```math
\tanh K_i=e^{-2ag_i}\le e^{-2ah},\qquad \tanh(aJ_{ij})\le aJ_{ij}.
```

These identities are exact at fixed a. The classical model therefore supplies the FKG
inequality for the time-dependent diagonal observables of $T_a$.

### 2. Summing paths after summing temporal segments

For a finite zero-field ferromagnetic Ising model, let $w_e=\tanh K_e$. Its high-temperature
expansion implies

```math
\langle\sigma_u\sigma_v\rangle
 \le\sum_{\gamma:u\to v\ \mathrm{simple}}\prod_{e\in\gamma}w_e.
```

Indeed an edge set with odd boundary $\{u,v\}$ contains such a path. Removing the path
leaves an even edge set disjoint from it. Dropping the disjointness restriction bounds
the remaining sum by the full even-edge partition sum, which cancels the denominator.
The inequality also holds for $u=v$ upon including the empty path.

Lift each path on the time cylinder to $G\times\mathbb Z$. Between consecutive spatial
edges its time coordinate is monotone: a reversal would immediately revisit a vertex.
Dropping all other self-avoidance restrictions gives an upper bound. A temporal segment
of displacement r has weight at most

```math
b_a(r)=e^{-2ah|r|},\qquad r\in\mathbb Z.
```

For paths with m spatial edges, summing their spatial weights is at most $(as)^m$.
Summing the intermediate time coordinates gives the convolution $b_a^{*(m+1)}$.
Thus on the infinite time cylinder the two-point bound is dominated by

```math
\sum_{m\ge0}(as)^m b_a^{*(m+1)}(n).
```

The passage from a finite time cylinder is made below; no infinite-volume assumption is
needed at this step. Fix $0<\kappa<2h$. The *signed* exponential weight is essential:

```math
B_a(\kappa)=\sum_{r\in\mathbb Z} e^{\kappa ar}b_a(r)
 =\frac{\sinh(2ah)}{\cosh(2ah)-\cosh(\kappa a)}.
```

Since $\sup_r e^{\kappa ar}b_a(r)\le1$, the convolution inequality between the weighted
supremum and sum norms gives

```math
b_a^{*(m+1)}(n)\le B_a(\kappa)^m e^{-\kappa a|n|}.
```

For negative n use the weight with the opposite sign; its sum is the same. If
$asB_a(\kappa)<1$, the geometric sum is therefore bounded by
$C_a e^{-\kappa a|n|}$, where $C_a=[1-asB_a(\kappa)]^{-1}$.
On a cylinder of circumference M, lifting the endpoint adds only the image sum

```math
\langle\sigma_{i,0}\sigma_{j,n}\rangle_M
 \le C_a\sum_{k\in\mathbb Z}e^{-\kappa a|n+kM|}.
```

For fixed n this tends to $C_a e^{-\kappa a|n|}$ as $M\to\infty$. The left side tends
to the ground-state correlation for $H_a$, by the Perron eigenvalue of $T_a$. Finally,

```math
\lim_{a\downarrow0}asB_a(\kappa)=\frac{4hs}{4h^2-\kappa^2}<1
 \quad\Longleftrightarrow\quad \kappa^2<4h(h-s).
```

Every $\kappa<2\sqrt{h(h-s)}$ therefore gives this decay bound for all sufficiently
small a. The allowed a depends only on h, s and $\kappa$, not on graph size.

### 3. From two-point correlations to the full spectral gap

It is necessary to control all diagonal observables, since a single spin channel need
not detect the first excitation. Here FKG supplies the missing implication directly.
For a real function f of one spin slice, put

```math
\ell_i(f)=\frac12\max_\sigma|f(\sigma^i)-f(\sigma)|,
\qquad L_f(\sigma)=\sum_i\ell_i(f)\sigma_i.
```

Both $L_f+f$ and $L_f-f$ are increasing functions. Apply FKG to the four pairs
$(L_f\pm f,L_q\pm q)$ on two time slices. Adding the two equal-sign inequalities and
then the two opposite-sign inequalities cancels the mixed covariances and gives

```math
|\mathrm{Cov}(f,q)|\le\mathrm{Cov}(L_f,L_q)
 =\sum_{i,j}\ell_i(f)\ell_j(q)\mathrm{Cov}(\sigma_i,\sigma_j).
```

The finite-cylinder spin means vanish by spin-flip symmetry. Taking $M\to\infty$
preserves this covariance inequality and the preceding two-point bound. In particular,

```math
0\le\langle f(0)f(na)\rangle_0-\langle f\rangle_0^2
 \le C_a\Big(\sum_i\ell_i(f)\Big)^2e^{-\kappa an}.
```

Let $\psi_0>0$ and $\psi_1$ be real normalized ground and first-excited eigenvectors of
$H_a$. The diagonal function $f(\sigma)=\psi_1(\sigma)/\psi_0(\sigma)$ is well defined
on the finite configuration space. It satisfies $f\psi_0=\psi_1$ and
$\langle f\rangle_0=0$. Its covariance is exactly $e^{-na\Delta(H_a)}$.
The last inequality, as $n\to\infty$, proves $\Delta(H_a)\ge\kappa$.
The prefactor may depend on the graph and a; this does not affect the decay rate.

Taking $a\downarrow0$ and using finite-dimensional eigenvalue continuity gives
$\Delta(H)\ge\kappa$. Taking $\kappa\uparrow2\sqrt{h(h-s)}$ proves the theorem.
Uniqueness also follows directly from the strictly positive entries of $e^{-tH}$ for
$t>0$, since every $g_i>0$ connects the configuration hypercube.

The only external correlation theorem used here is
[Fortuin, Kasteleyn and Ginibre (1971)](https://doi.org/10.1007/BF01651330).
The path estimate, temporal sum, covariance reduction, and spectral argument are given
above. This graph theorem is a mathematical proof; it has not been formalized in Lean.
The numerical checker verifies the exact rational evaluations of its explicit bound.
It does not check the FKG theorem or the proof in this document.

## Verification scope

The fast checker verifies model-local source hashes, exact row models and RP factors, stored
spectrum hashes and dimensions, fresh interval partition sums, and the exact cube chain.
The full checker also rebuilds every symmetry block and verifies its complete spectrum.
The spectral trusted base is IEEE binary64 round-to-nearest with gradual underflow, the
standard dot-product rounding bound for BLAS (including arbitrary summation order and FMA),
and interval arithmetic. Unverified eigenvalues never enter the certified defects.

The transfer argument and the finite asymmetric recurrence are mathematical proofs supplied
here and in `paper/thermal-mpo-gap-2plus1d.tex`. The general small-defect theorem has the Lean
formalization stated above. The complete Ising certificate is checked by Julia; there is
no claim of a Lean proof from the Ising Hamiltonian through its numerical seeds. The
bootstrap certificate covers periodic square tori. The graph theorem additionally covers
open boundaries and arbitrary finite graphs satisfying its incident-coupling condition.
Neither argument determines the critical field.

## Recorded verification

Julia 1.12.4, with one BLAS thread, produced complete spectra of dimensions 16, 256, and
65536. These are represented by 8, 12, and 24 blocks, respectively, with exact multiplicities.
The largest block has dimension 2392. At fields 18, 20, and 24 the largest eigenvalue radii
are below $4.21\times10^{-7}$, $4.70\times10^{-7}$, and $5.68\times10^{-7}$, respectively.
All three complete spectral replays imply the stored seed inequalities, with uniform gap
lower bounds $0.5943744623$, $2.1736330392$, and $4.8055785105$ (rounded down).

The 42 model tests include an independent bit-basis construction of the small tori, the
identity $\mathrm{Tr}\,H_{4,4}^2=2^{16}(16g^2+32)$ checked against interval spectral
moments, exact recurrence and tail checks, and deliberate corruption of the advertised
gap, seed defects, and spectral hashes. An unchanged temporary copy is checked before its
altered versions, so these rejection tests do not depend on a missing environment.

The core's complete `Pkg.test()` suite passed after the sparse-permutation change. That
change replaces only the storage of exact signed permutation matrices in symmetry and
reflection checks; it changes neither the spectral verification inequalities nor the
arithmetic model. Existing model checkers and tests were also rerun as a regression check.
Lossless compression of the stored records was separately checked to preserve every seed,
recurrence step, and spectral record.

The refinement adds sixteen tests, bringing the total to 58. All passed, including rejection
of an altered refined gap and an altered ground-energy enclosure. Regeneration preserved
the original seed, chain, and spectral records by exact equality for all three fields.
Full spectral replay verified both sets of seed bounds in 217.1 seconds over all three fields.
`lake build CertifiedGaps` passed, and the three new principal declarations
(`cube_growth`, `cascade`, and `gapGT_of_cube_seeds`) depend only on `propext`,
`Classical.choice`, and `Quot.sound`.

For the graph theorem, 25 arithmetic and numerical checks passed, including independent
Hamiltonians on eighteen weighted graphs with nonuniform transverse fields. These numerical
tests are checks on the implementation, not part of the proof. The existing 58 bootstrap
checks also passed. After adding the new sources, all original seed, chain, refinement,
spectral, and finite-volume records were preserved by exact equality; their metadata was
regenerated and the fast checkers passed. No spectral computation is used to evaluate the
analytic bound. The separate analytic checker completed in 2.35 seconds including startup.
