# A finite certificate for the square-lattice transverse-field Ising model

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

The uniform bounds remain well below the sixteen-spin gaps. At $g/J=24$, for example,
the finite gap exceeds $43.91119J$, while the uniform bound is $5.52173J$. Spectral enclosure
errors are below $10^{-6}J$ per eigenvalue. Sharper eigenvalue arithmetic alone would thus
have little effect. Larger spatial seeds or sharper directional estimates are needed for
a substantial further increase.

## The ordered side of the transition

For the Pauli normalization used here, numerical studies place the square-lattice
transition at $g/J\simeq3.04438$; see
[Schmitt et al., Science Advances 8, eabl6850 (2022)](https://pmc.ncbi.nlm.nih.gov/articles/PMC9481121/).
The appropriate low-field target on symmetric periodic tori is an isolated pair of
low-energy states. With eigenvalues counted with multiplicity, define

```math
\delta_L=E_1(L)-E_0(L),\qquad \gamma_L=E_2(L)-E_1(L).
```

A certificate for this phase should control both the vanishing splitting $\delta_L$ and a
positive uniform lower bound on $\gamma_L$. The existing `GapGT` conclusion concerns
$\delta_L$ and therefore cannot establish this different statement.

At $g=0$, the two constant spin configurations have energy $-2JL^2$. Every nonconstant
configuration has at least four disagreeing edges on a square torus with $L\ge4$, and each
such edge costs $2J$. A single flipped spin attains this minimum. Consequently
$\delta_L=0$ and $\gamma_L=8J$ exactly.

There is an exact obstruction to the present seed criterion at this endpoint. Global spin
reversal pairs every classical configuration with a distinct configuration of the same
energy. Summing one representative per pair gives

```math
Z_L(\beta)=2\sum_\alpha e^{-\beta e_\alpha},\qquad
P_\tau(L,\beta)=\frac12
 \frac{\sum_\alpha e^{-2\beta e_\alpha}}
 {(\sum_\alpha e^{-\beta e_\alpha})^2}\le\frac12.
```

Thus the thermal defect is at least $1/2$ at every temperature. At positive field, when
higher levels have negligible thermal weight, the two-level approximation instead gives

```math
P_\tau(L,\beta)\simeq
\frac{1+e^{-2\beta\delta_L}}{(1+e^{-\beta\delta_L})^2}.
```

This is close to $1/2$ when $\beta\delta_L\ll1$. For fixed finite volume and nonzero
splitting, it tends to one only when the temperature resolves the splitting. A uniform
gap above the pair cannot be inferred from such a resolution scale.

Complete floating-point spectra on the 4 by 4 torus give:

| $g/J$ | $E_0/(16J)$ | $\delta_4/J$ | $\gamma_4/J$ |
|---:|---:|---:|---:|
| 1 | -2.12566234693 | $5.51378\times10^{-7}$ | 7.28553502 |
| 2 | -2.51188715235 | 0.0150166052 | 4.71774367 |

The first four levels, in units of $J$, are approximately
$(-34.01059755084623,-34.01059699946832,-26.72506197570735,-26.72419743700171)$
at $g/J=1$, and
$(-40.19019443767035,-40.17517783245594,-35.45743416578414,-34.50218873303338)$
at $g/J=2$. These are unverified LAPACK results, with no eigenvalue enclosures or
thermodynamic extrapolation. Reproduce them from the repository root with:

```sh
OPENBLAS_NUM_THREADS=1 julia --project=code/julia -e '
using CertifiedGaps.VerifiedSpectra: float_ring_spectrum
include("code/models/square_tfim/src/SquareIsing.jl")
using .SquareIsing
for g in (1, 2)
    seconds = @elapsed E = float_ring_spectrum(row_model(4, g), 4).E
    println((g=g, levels=E[1:4], density=E[1]/16, splitting=E[2]-E[1],
        gap_above_doublet=E[3]-E[2], seconds=seconds))
end
'
```

The eigensolves took 27.4 and 25.9 seconds respectively on the development machine with
one BLAS thread, excluding initial package startup.

A possible extension would use the two spin-flip sectors separately. Their thermal traces
are positive, but spatial doubling of globally projected traces requires a new gluing
argument. Another possibility is to control an isolated two-dimensional low-energy space
directly. Purity near $1/2$ alone does not certify this: the three probabilities
$(2/3,1/6,1/6)$ have purity exactly $1/2$. In particular, substituting $2P_\tau$ into the
unique-state recurrence has no justification. No ordered-phase extension is proved here.

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
of `WindowInput` from spatial transfer moments and the numerical seed arithmetic remain
outside Lean. The finite asymmetric prefix is covered by the following theorem.

## General finite-chain criterion

Set $Q(a,b,c)=(1-a)^4(1-b)^2(1-c)$ and $f(z)=z^2/(2(1-z)^2)$.
Let $(a_j,b_j,c_j)$, $0\le j\le K$, be upper bounds on the three staggered
seed defects. Suppose each bound is at most one, and each step satisfies

```math
\begin{aligned}
a_{j+1}&\ge f(1-Q(a_j,b_j,c_j)),\\
b_{j+1}&\ge f(1-Q(b_j,c_j,a_{j+1})),\\
c_{j+1}&\ge f(1-Q(c_j,a_{j+1},b_{j+1})).
\end{aligned}
```

Each of the three products appearing here must be positive. This permits rational upper
rounding at every update. Suppose the final three defects are bounded by
$0<\varepsilon\le10^{-3}$. If $r>1$ satisfies

```math
r^{2^j}(1-Q(a_j,b_j,c_j))<Q(a_j,b_j,c_j)\quad(0\le j<K),
\qquad r^{2^K}25\varepsilon\le1,
```

then every even square side $L\ge n_0$ has a unique ground state and
$\Delta_L>\log(r)/\beta_0$, under the scalar moment and window assumptions above.

**Proof.** The lower bound $Q(a_j,b_j,c_j)$ on the cube product bounds the last
purity along any reordered path through that cube. Monotonicity of $f$ and the
single-direction squaring inequality give the first update. Two cyclic permutations give
the other updates. Induction proves the three purity bounds through step $K$.
For a side length in an earlier dyadic window, the displayed strict inequality gives the
gap directly from thermal purity. For $L\ge n_0 2^K$, shift the data to scale $K$ and
apply the small-defect theorem. The tail inequality ensures
$\log(r)/\beta_0\le\log(1/(25\varepsilon))/(2^K\beta_0)$.

This is `Spacetime.gapGT_of_cube_chain` in
[`CubeChain.lean`](../../lean/CertifiedGaps/CubeChain.lean). The Lean statement allows
$r>0$; $r>1$ is the case giving a positive gap. Its hypotheses on the numerical chain
are finite rational inequalities when the supplied data are rational. No Ising-specific
high-field theorem is used. The supplied Julia records have not been converted into
Lean witnesses, so the numerical certificates still rely on the Julia checker.

## Verification scope

The fast checker verifies model-local source hashes, exact row models and RP factors, stored
spectrum hashes and dimensions, fresh interval partition sums, and the exact cube chain.
The full checker also rebuilds every symmetry block and verifies its complete spectrum.
The spectral trusted base is IEEE binary64 round-to-nearest with gradual underflow, the
standard dot-product rounding bound for BLAS (including arbitrary summation order and FMA),
and interval arithmetic. Unverified eigenvalues never enter the certified defects.

The transfer argument is supplied here and in `paper/thermal-mpo-gap-2plus1d.tex`.
The general finite-chain criterion has the Lean formalization stated above. The complete
Ising certificate is checked by Julia; there is
no claim of a Lean proof from the Ising Hamiltonian through its numerical seeds. The
certificate covers periodic square tori, with no conclusion about the critical field.

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

The finite-chain extension adds `cube_step`, `cube_purities_ge_chain`, and
`gapGT_of_cube_chain`. A full `lake build CertifiedGaps` passed. An explicit axiom audit of
these declarations reports only `propext`, `Classical.choice`, and `Quot.sound`.
After removing the separate high-field calculation from the default drivers, all three
bootstrap certificates were regenerated and checked. Their original seeds, chains, refined
seeds and chains, spectral records, and finite-volume enclosures were preserved by exact
equality. Source metadata changed to describe the active bootstrap computation.
The 58 model regression tests passed again, including rejection of altered certificate
claims; the test sets took 130.2 seconds, excluding startup.
