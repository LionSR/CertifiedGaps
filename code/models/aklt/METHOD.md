# On-site reflection obstruction and conditional trajectory bounds

This is a standalone research note, not an AKLT gap certificate or an amendment to the
initial manuscript. It proves an obstruction to every on-site antiunitary reflection for
the spin-3/2 AKLT bond, then records elementary error-tolerant bootstrap criteria. **The
model-specific thermal gluing estimate needed to apply those criteria is still unproved.**

The discussion refines the existing [2+1D companion note][companion].
It does not change that note, any certified bound, or the published paper.

## 1. Conventions

The physical single-site space is $`\mathcal H=\mathbb C^4`$, carrying spin $`3/2`$.
The bond interaction $`h=P^{(3)}`$ is the orthogonal projector onto total spin three.
All traces below are ordinary, unnormalized traces unless a normalized trace is explicitly
written. The Hilbert-Schmidt norm is $`\|A\|_2^2=\mathrm{Tr}(A^\dagger A)`$.

For any antiunitary $`J:\mathcal H\to\mathcal H`$, set $`\theta_J(A)=JAJ^{-1}`$.
The necessary reflected-positivity test considered here is

```math
\mathrm{Tr}\!\left[e^{-t h}(A\otimes\theta_J(A))\right]\ge0
\qquad\text{for every }A\in\mathcal B(\mathcal H).
```

A negative value for even one Hermitian test matrix disproves this condition. The theorem
below allows arbitrary $`J`$, not only those compatible with an involutive full reflection.
It therefore includes, but is not restricted to, standard on-site reflection constructions.

## 2. An antiunitary fixed-point lemma

**Lemma 1.** Every antiunitary on a finite-dimensional complex Hilbert space of dimension
$`d\ge3`$ fixes a nonzero traceless Hermitian operator under conjugation.

**Proof.** The linear operator $`J^2`$ is unitary. Choose an eigenvector $`v\ne0`$ with
$`J^2v=\lambda v`$. The complex space $`V=\mathrm{span}_{\mathbb C}\{v,Jv\}`$ has dimension
one or two, and $`J(V)=V`$. An antiunitary also preserves orthogonality, so the orthogonal
projector $`P_V`$ is fixed by conjugation with $`J`$. Since $`0<r=\dim V<d`$,
$`A=P_V-(r/d)I`$ is nonzero, Hermitian, traceless, and fixed. This argument does not assume
$`J^2=\pm I`$. ∎

## 3. The coherent-state frame is strictly positive

Let $`\rho_\Omega=|\Omega\rangle\langle\Omega|`$ be the spin-3/2 coherent-state projector,
and normalize spherical measure by $`\int d\Omega=1`$. Parallel coherent vectors on the
two sites lie in the total-spin-three representation. Its coherent-state resolution gives

```math
P^{(3)}=7\int\rho_\Omega\otimes\rho_\Omega\,d\Omega.
```

Consequently, on Hermitian matrices,

```math
B(A,A):=\mathrm{Tr}[P^{(3)}(A\otimes A)]
       =7\int[\mathrm{Tr}(\rho_\Omega A)]^2\,d\Omega.
```

This form is strictly positive when $`A\ne0`$. For completeness, in a stereographic chart
an unnormalized coherent vector has coordinates
$`(1,\sqrt3 z,\sqrt3 z^2,z^3)`$. Its expectation of $`A`$ is a polynomial with coefficients
$`\sqrt{\binom3a\binom3b}\,A_{ab}`$ multiplying $`\bar z^a z^b`$. These monomials are
linearly independent, because $`z=x+iy`$ and $`\bar z=x-iy`$ give an invertible linear
change of polynomial variables. A zero symbol therefore implies $`A=0`$.

**Lemma 2.** The frame operator

```math
K(A)=7\int\rho_\Omega\,\mathrm{Tr}(\rho_\Omega A)\,d\Omega
```

has eigenvalues $`7/4,21/20,7/20,1/20`$ with multiplicities $`1,3,5,7`$. In particular,
$`B(A,A)\ge\|A\|_2^2/20`$ for every Hermitian $`A`$.

**Proof.** The coherent-symbol map and its adjoint identify the nonzero frame spectrum with
that of the spherical integral kernel $`7|\langle\Omega|\Omega'\rangle|^2`$. Its overlap is
$`|\langle\Omega|\Omega'\rangle|^2=[(1+\Omega\cdot\Omega')/2]^3`$. The spherical-harmonic
decomposition, equivalently the Funk-Hecke formula with normalized measure, gives

```math
k_\ell=\frac72\int_{-1}^1\left(\frac{1+x}{2}\right)^3P_\ell(x)\,dx,
\qquad \ell=0,1,2,3.
```

Inserting the four Legendre polynomials yields the displayed rational numbers; their
multiplicities are $`2\ell+1`$ and sum to sixteen. The frame is self-adjoint in the physical
Hilbert-Schmidt inner product. The smallest eigenvalue proves the inequality. ∎

Here $`K`$ is a coherent-state frame operator. It must not be confused with a Choi matrix
or a differently transposed realignment of the bond. The exact checker independently
verifies the projector, the frame spectrum, and the physical inner product in its gauge.

## 4. No on-site antiunitary makes the bond reflection positive

**Theorem 3.** For every antiunitary $`J`$ on $`\mathbb C^4`$, there is a traceless
Hermitian $`A`$ of Hilbert-Schmidt norm one, independent of $`t`$, such that for all $`t>0`$,

```math
\mathrm{Tr}\!\left[e^{-tP^{(3)}}(A\otimes\theta_J(A))\right]
\le-\frac{1-e^{-t}}{20}<0.
```

**Proof.** Normalize the fixed Hermitian operator from Lemma 1. Since
$`e^{-tP^{(3)}}=I-(1-e^{-t})P^{(3)}`$ and $`\mathrm{Tr}A=0`$, the left side equals
$`-(1-e^{-t})B(A,A)`$. Apply Lemma 2. ∎

Thus choosing a different on-site antiunitary or merely changing the local physical bases
does not rescue this bond construction. The older quadrupole witness
$`Q=\mathrm{diag}(1,-1,-1,1)`$ has $`B(Q,Q)=7/5`$ and is the special case for physical
spin time reversal. The qualitative proof extends to the maximal-total-spin projector
for every spin $`s\ge1`$, by the same coherent-state resolution and Lemma 1. Spin one-half
is excluded: a two-dimensional Kramers antiunitary can lack a proper invariant subspace.

### A finite-lattice consequence with explicit geometry

**Corollary 4.** Let $`H=\sum_{e\in E}P^{(3)}_e`$ on a finite loopless multigraph with
$`N`$ sites and $`M=|E|\ge1`$ edges, counted with multiplicity. Suppose a reflection exchanges
two disjoint halves by on-site antiunitaries and exchanges distinct sites $`u,v`$ connected
by $`k\ge1`$ edges. Then its thermal functional fails reflection positivity whenever

```math
0<\beta\le\frac1{80M^2}.
```

**Proof.** Use the normalized fixed witness from Theorem 3 for the antiunitary identifying
the two sites, with identity elsewhere. Write $`W=A_u\theta_J(A)_v`$ and
$`\tau_N(X)=4^{-N}\mathrm{Tr}X`$. Then

```math
f(0)=0,\qquad
f'(0)=-\tau_N(HW)=-\frac{k B(A,A)}{16}\le-\frac1{320},
\qquad f(\beta)=\tau_N(e^{-\beta H}W).
```

An edge other than the $`k`$ edges joining $`u,v`$ leaves at least one traceless witness
unpaired, so its trace vanishes. Since $`0\le H\le MI`$, scalar functional calculus gives

```math
0\le e^{-\beta H}-I+\beta H\le\frac{\beta^2M^2}{2}I.
```

Also $`\tau_N(|W|)=(\mathrm{Tr}|A|)^2/16\le1/4`$ by Hilbert-Schmidt Cauchy-Schwarz.
The remainder in $`f`$ is therefore at most $`\beta^2M^2/8\le\beta/640`$, giving
$`f(\beta)\le-\beta/640<0`$. Dividing by the positive partition function does not change
the sign. ∎

The graph and interaction assumptions matter. This statement is not asserted for arbitrary
multisite cross-plane terms or reflections acting nonlocally inside a block. It obstructs
an all-temperature RP construction, not RP restricted to a low-temperature set.

## 5. A conditional absolute-error envelope

For frustration-free finite Hamiltonians with at least one exact zero-energy state, set
$`Z(n,m,\beta)=\mathrm{Tr}\,e^{-\beta H_{n,m}}`$ and $`R=Z-1`$. Along the dyadic sequence
$`(n_j,m_j,\beta_j)=2^j(n_0,m_0,\beta_0)`$, write $`R_j=R(n_j,m_j,\beta_j)`$. Define

```math
\eta_j=\max\{0,\log Z(2n_j,2m_j,2\beta_j)
                    -4\log Z(n_j,m_j,2\beta_j)\}.
```

Temporal squaring, obtained by squaring the nonnegative excited Boltzmann weights, gives
$`R(n_j,m_j,2\beta_j)\le R_j^2`$. Thus, without any RP assumption,

```math
R_{j+1}\le e^{\eta_j}(1+R_j^2)^4-1.
```

**Theorem 5.** Let $`0<q<1`$. If a seed and a uniform model-specific error estimate give

```math
R_0\le\frac q{10},\qquad
\eta_j\le\frac1{40}q^{2^{j+1}}\quad\text{for every }j\ge0,
```

then $`R_j\le q^{2^j}/10`$, the ground state is unique at every dyadic size, and
$`\Delta_{n_j,m_j}>\log(1/q)/\beta_0`$.

**Proof.** For $`0\le r\le1/10`$, $`0\le\eta\le1/40`$, use
$`e^\eta\le1/(1-\eta)\le40/39`$ to obtain

```math
e^\eta(1+r^2)^4-1\le5r^2+2\eta.
```

Indeed, the coefficient of $`r^2`$ is at most
$`(40/39)(4+6/100+4/10000+1/1000000)=1353467/325000<5`$; the coefficient of $`\eta`$
is at most $`40/39<2`$. Inserting the asserted envelope yields exactly

```math
5\left(\frac{q^{2^j}}{10}\right)^2
 +2\frac{q^{2^{j+1}}}{40}
 =\frac{q^{2^{j+1}}}{10}.
```

This proves the induction. Since $`R_j<1`$, additional zero eigenvalues are impossible.
Finally $`e^{-\beta_j\Delta_j}\le R_j\le q^{2^j}/10`$, so
$`\Delta_j\ge\log(1/q)/\beta_0+\log10/\beta_j`$. ∎

The hypothesis uses an absolute envelope, not an inequality relative to the unknown actual
$`R_j^2`$. It needs only the specified tuples $`(n_j,m_j,\beta_j)`$,
$`(n_j,m_j,2\beta_j)`$, and $`(2n_j,2m_j,2\beta_j)`$. In particular, there is a doubled-temperature
companion sequence; this is not a claim to use only one ray in size-temperature space.
No all-temperature inequality at each fixed volume is required. No interpolation to
intermediate periods follows without a separate hypothesis.

**Corollary 6.** Suppose $`0<R_0<1/10`$ and a proved estimate gives
$`\eta_j\le a e^{-\mu(2^j-1)}`$ with $`0<a<1/40`$, $`\mu>0`$. Any positive

```math
\kappa\le\min\left\{\log\frac1{10R_0},\frac\mu2,
                       \frac12\log\frac1{40a}\right\}
```

gives $`\Delta_j>\kappa/\beta_0`$.

**Proof.** Set $`q=e^{-\kappa}`$. The seed condition is immediate, and
$`\eta_j e^{2\kappa2^j}\le a e^{\mu-(\mu-2\kappa)2^j}\le a e^{2\kappa}\le1/40`$.
Apply Theorem 5. ∎

Thus spatial-side-length exponential errors suffice, with their initial constants controlled.
This is a conditional criterion, not a claim that AKLT thermal gluing errors obey it.

## 6. Scale-dependent doubling coefficients

**Proposition 7.** Suppose $`R_{j+1}\le C_jR_j^2`$ with $`C_j>0`$, and define

```math
S=\sum_{j=0}^\infty2^{-j-1}\log\max\{1,C_j\}.
```

If $`S<\log(1/R_0)`$, then

```math
R_j\le\exp\{-2^j[\log(1/R_0)-S]\},\qquad
\Delta_j\ge\frac{\log(1/R_0)-S}{\beta_0}.
```

**Proof.** Take logarithms of the recurrence, divide by $`2^{j+1}`$, and telescope.
The partial weighted sum is at most $`S`$. The resulting bound is less than one, so uniqueness
and the gap estimate follow as above. ∎

For example, $`C_j\le C_*2^{pj}`$ with $`C_*\ge1`$, $`p\ge0`$ costs only
$`S\le\log C_*+p\log2`$. Unlike the gluing-error formulation, the actual coefficient
$`C_j=R_{j+1}/R_j^2`$ involves only the main dyadic sequence. Its model-specific control
remains necessary; finitely many measured coefficients prove no uniform statement.

## 7. Exact verification and its limits

Run from the repository root, with Julia 1.12:

```sh
julia --startup-file=no code/models/aklt/scripts/verify_rp_trajectory.jl
```

This uses only Julia standard libraries and `Rational{BigInt}`, with no floating-point
eigensolver, external package, or downloaded evidence. The physical diagonal similarity
$`D=\mathrm{diag}(1,\sqrt3,\sqrt3,1)`$ makes the bond rational:

```math
\widetilde P_{(a,b),(c,d)}
 =\mathbf1_{a+b=c+d}\frac{\binom3c\binom3d}{\binom6{a+b}},
\qquad 0\le a,b,c,d\le3.
```

The checker verifies idempotence, rank, and an independent total-spin polynomial construction.
It builds the frame in matrix-unit coordinates and checks its positive weighted
Hilbert-Schmidt metric, exact spectral projections and multiplicities, the quadrupole trace,
a finite-lattice cancellation example, and rational constants and recurrence identities.
Deliberately corrupted matrices fail the tested identities. The same checks are included in
the existing AKLT test suite.

The universal quantifiers in Theorems 3 and 5 are established by the proofs above, not by
sampling antiunitaries, temperatures, or lattice sizes. Finite exact checks support the
algebra and constants; this is not a Lean formalization. The existing 8-spin thermal seed
and compressed thermal-MPO results remain uncertified floating-point diagnostics.

## 8. Prior art and novelty status

Targeted literature check, November 17, 2026:

- Jaffe and Janssens, [*Characterization of Reflection Positivity: Majoranas and Spins*][jj16],
  especially Theorems 1.1 and 5.2 and Section 6,
  characterize all-temperature RP through cross-plane couplings and analyze gauge changes.
  The small-temperature derivative obstruction is part of that established framework.
- Jaffe and Janssens, [*Reflection Positive Doubles*](https://arxiv.org/abs/1607.07126),
  especially Theorem V.10 and Section VII.2, give the broader algebraic characterization,
  including matrix-algebra reflections. Their scope is broader than antiunitary conjugation;
  the present theorem must not be described as excluding all such algebraic reflections.
- Amiet and Weigert, [*Discrete Q- and P-symbols for spin s*][aw00],
  establish coherent-projector informational
  completeness. That ingredient is standard, not a new contribution here.
- The repository's earlier [2+1D note](../../../paper/thermal-mpo-gap-2plus1d.tex) already
  contains the physical-time-reversal quadrupole obstruction, the exact perturbed excess
  update, and the uniformly bounded-coefficient criterion. Those are not claimed as new.

The independent derivations here extend the on-site witness to every antiunitary and give
explicit sufficient error envelopes and weighted scale-loss budgets. The targeted search
did not locate these exact combined statements in the checked sources. **This is not an
exhaustive novelty determination or a priority claim.** The scalar recurrence refinements
are elementary auxiliary estimates, not a substitute for the missing AKLT gluing theorem.

## 9. Remaining target

Prove a bound on the positive part of the combined thermal gluing defect at the specified
tuples, or prove a summable weighted-log bound for the actual doubling coefficients. The
ground-state PEPS boundary Gram operators and the fixed-width/depth projection filter are
different objects; neither presently supplies these thermal estimates. The results above
do not disprove scalar gluing, exclude blocked/nonlocal reflection, establish an AKLT gap,
or cover intermediate periods.

[jj16]: https://arxiv.org/abs/1506.04197
[aw00]: https://arxiv.org/abs/quant-ph/9906099
[companion]: ../../../paper/thermal-mpo-gap-2plus1d.tex
