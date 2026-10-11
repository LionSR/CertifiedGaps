# Independent mathematical and implementation audit

Reviewer: independent delegated review agent, 2026-10-11. This review did not
author the original calculation and did not run a duplicate twelve-spin job.
It inspected `rp.jl`, `common.jl`, `verified_blocks.jl`, `certify.jl`, both test
files, the stored delta=7/10 certificate, and the transfer, coupling,
interpolation and gap-recovery arguments in the repository notes. The numerical
review concerns the original certificate algorithm; the separate publication
workflow is adapting the core package name and replaying the calculation.

**Finding:** no substantive mathematical obstruction was found to the odd
three-cell seed. The positive transfer proof below supplies the extra hypothesis
that the generic self-adjoint-transfer theorem does not supply. The enclosure
algorithm has the required complete-spectrum and outward-bound structure.
Publication as a replay-verified result remains conditional on completion and
inspection of the fresh full spectral replay. This document alone does not
attest that replay has completed.

## 1. Exact model and transfer orientation

Write J-=1-delta and J+=1+delta. The integer Cartesian matrices in `rp.jl`
are real antisymmetric, with Sa=-i Ga, the spin-one commutation relation and
sum Sa^2=2I. Thus -J Sa tensor Sa=J Ga tensor Ga, with positive rational J
for every requested parameter. This is an exact identity, not a tolerance test.
For a cell with weak internal bond, D=-J- sum Ga tensor Ga is real symmetric.
The internal site swap P maps A=I tensor Ga to B=Ga tensor I and fixes D.
Consequently -h_inter=J+ sum A tensor (P conjugate(A) P).
The naive claim that the intercell operators are identical on both cell factors
would be false; the implemented internal-reflection identity has the correct
orientation.

Here is a direct transfer construction that avoids relying on that claim.
For J>0 let mu_J be the **unnormalised** measure on finite histories of labelled
events in [0,beta], with weight J^m on each ordered m-event simplex, and a unit
atom at the empty history. There are three labels. Its mass is exp(3 beta J).
For two histories a,u, define F(a,u) as the trace of the time-ordered product
of the Ga insertions from their union. Coincident event times have measure zero
and may be assigned a consistent order. Since both incident bonds insert exactly
the same Ga on this site, exchanging their history names changes no ordered
product: F(a,u)=F(u,a). All Ga are real, so F is real. Since ||Ga||=1,
|F|<=3, including the empty product.

Set H+=L2(mu_J+) and H-=L2(mu_J-). The integral operator K:H- -> H+ with
kernel F is Hilbert-Schmidt since both measures are finite and F is bounded.
Integrating the weak history inside a cell gives the strong-boundary kernel

    X(a,b) = integral F(a,u) F(u,b) dmu_J-(u) = (K K*)(a,b).

Therefore X is positive and trace class. Expanding the bounded finite-ring
Hamiltonian exponential into bond events and taking the physical trace gives
exactly Z_L(beta)=Tr X^L. Independent bond histories are equivalent to the
global ordered Dyson expansion; there is no Poisson normalisation factor
because the measures above are unnormalised. Absolute convergence follows from
the bounded insertions and finite history measures. This identity holds for
every positive integer number of cells (with both periodic bonds retained for
the one-cell ring); only L>=3 is needed here. In particular all transfer
eigenvalues lambda_i(beta) are nonnegative at every beta>0.

## 2. Why n0=3 has the same recurrence and uniform conclusion

For every integer n>=1, X^n/Z_n is a density operator. Its purity is
P_s(n,beta), so P_s<=1 and the spatial squaring lemma applies at odd n.
Thermal density operators give P_t<=1 without any parity assumption.
Partition-function cancellation gives the two coupling identities in the note
for every such n. Define

    z(s,t)=1-(1-s)^2(1-t),    f(z)=z^2/[2(1-z)^2].

Starting from defects bounded by (s,t), the first coupling gives a defect at
(n,2 beta) bounded by z(s,t). Spatial squaring gives
s'=f(z(s,t)) at (2n,2 beta). The second coupling gives a thermal defect at
(4n,beta) bounded by z(t,s'), and temporal squaring gives
t'=f(z(t,s')). This is precisely the unequal-defect recurrence used by the
certified implementation. Its derivation has no remaining parity restriction.

For n_j=2^j n0, beta_j=2^j beta0, and n_j<=L<=2n_j, nonnegative transfer
eigenvalues and monotonicity of sequence p-norms imply

    Z_L(beta_j) <= Z_nj(beta_j)^(L/n_j),
    Z_L(2 beta_j) >= Z_2nj(2 beta_j)^(L/(2 n_j)).

Hence, writing alpha=L/(2n_j)<=1,

    P_t(L,beta_j)
      >= [P_t(2n_j,beta_j) P_s(n_j,beta_j)^2]^alpha
      >= (1-t_j)(1-s_j)^2 = 1-z(s_j,t_j).

The certified orbit requires z(s_j,t_j)<1/2 at every scale. Since the largest
Gibbs eigenvalue w0 is at least its purity, w0>1/2. A degenerate ground state
would have w0<=1/2, proving uniqueness without an external uniqueness
assumption. Moreover

    exp(-beta_j Delta_L) <= (1-w0)/w0
                        <= z(s_j,t_j)/(1-z(s_j,t_j)).

Thus Delta_L>=b_j, where b_j is the rate bounded by the core recurrence's
admissibility and tail theorem. Its certified lower bound on inf_j b_j works
for every L>=n0 by selecting its enclosing dyadic interval. Restricting this
conclusion to the requested even cell lengths gives L=4,6,8,... for n0=3,
or N=8,12,16,... physical spins. It does not accidentally claim N=6 or all
even physical spin counts. No extension of the generic theorem has been
assumed without supplying the additional positivity and interpolation argument.

## 3. Complete spectra and numerical verification

The Sz Hamiltonian uses integer digits m=-1,0,1. Each allowed nearest-neighbour
exchange has coefficient J: the spin-one square roots cancel the 1/2 in the
ladder expression. Diagonal entries are the exact weighted m_j m_(j+1) sums.
For the specified sizes and rational parameters the integer entries fit exactly
in Float64. Charge reversal, two-site translation, and reversal of site order
preserve this Hamiltonian. Reversal preserves the alternating bond parity.

The character projectors use exact integer twice-cosines at cell periods
2,3,4,6; symmetry, idempotence after scaling, and trace/rank are asserted per
orbit. Positive and negative generic momenta have equivalent spectra. Restricting
to reflection + takes one copy, with multiplicity two: on a generic momentum
pair, T-T^(-1) gives an invertible intertwiner between the two reflection
parities. Both reflection parities are retained at k=0 and k=L/2. Spin flip
supplies multiplicity two for M>0 and is split explicitly at M=0. These are
exact representation facts, not numerical near-degeneracy assumptions.

An independent JSON read found 7 small-ring Sz blocks with weighted dimension
729 and maximum dimension 141, and 79 large-ring blocks with weighted dimension
531441 and maximum dimension 11596. These agree with 9^3 and 9^6. Dimension
counts alone would not prove completeness, but together with these disjoint
charge/character projectors and the exact rank assertions they do.

The polar/Weyl argument in the report is valid. If ||Q'Q-I||<=e<1 and
||HQ-Q Lambda||<=r, then Q=US has invertible S, sigma_min(S)>=sqrt(1-e),
and ||S-I||<=e/(1+sqrt(1-e)). Comparing H with U Lambda U' gives the stated
radius (r+2||Lambda|| ||S-I||)/sqrt(1-e). Weyl bounds every ordered eigenvalue
with multiplicity, including clusters. It is not merely a residual test for
isolated eigenvectors.

For symmetry bases, the independently checked projector residual controls PB-B;
the Gram defect of PB is at most e+2||B||q+q^2. The polar basis of PB spans the
exact projector range because its dimension equals the exact rank and this
defect is below one. Disjoint orbits make the full operator-norm distance the
maximum orbit distance. The bound 4L on ||H|| is valid: ||Sa.Sa||=2 for a
spin-one bond and the sum of alternating positive couplings is 2L.
The projected-matrix construction allowance and eigenvector allowance therefore
refer to the exact model, rather than silently treating numerical bases as exact.

Residuals and Gram defects are accumulated in Float64 and enlarged with normwise
IEEE error bounds; this is the allowed rigorous-norm-bound route, not high
precision matrix multiplication. The gamma factors conservatively include dot
products, multiplication, subtraction and scaling. Norm sums include an absolute
underflow allowance. Scalar bounds use directed MPFR or interval arithmetic,
including lower-rounded polar denominators. The projection bounds use substantial
gamma-factor slack. This audit assumes standard IEEE arithmetic, no fast-math
reassociation invalidating the error model, and the stated single-thread BLAS
configuration; it is not a formal verification of BLAS or MPFR. No substantive
direction error was found in these bounds.

## 4. Partition directions, stored seed, and audit limits

Each exact spectrum lies within the computed centres plus/minus rational radii.
Exponentials, summation, logarithms, and the common energy shift are enclosed by
IntervalArithmetic. For the spatial purity the safe exponent is lower(logZ_2n)
minus twice upper(logZ_n); for thermal purity it is lower(logZ_2n_2beta) minus
twice upper(logZ_2n_beta). The implementation uses these correct directions
and rounds the resulting defect upper endpoint upward to a rational.

The independently read original stored inputs are

    epsilon_s = 58597799297//500000000000
    epsilon_t = 48363543279//1000000000000
    beta0     = 17//4

with stored bound 0.2187941698399417. The decimal lower claim
0.218794169839 is weaker and correctly rounded downward. Numerical reference
gaps, scan rankings, and extrapolated seeds remain numerical observations.

The fast checker verifies exact RP, rational-seed calculations, counts and source
provenance, plus small numerical comparisons. It does not reconstruct twelve-spin
eigenvalue radii from model inputs. Source hashes attest provenance, not correctness;
the stored `positive_transfer_gram_proof` boolean is metadata, not a proof checker.
The mathematical argument above supplies that proof. No Lean formalisation of this
model, transfer construction, or numerical enclosure algorithm is claimed.

A fresh complete `certify.jl` replay is required to close the independent
computational validation requested for publication. Preserve the original JSON
and identify the replay result, runtime, source hashes and any seed/bound drift
separately. Updating package-import compatibility must not silently replace the
historical sources against which the original results were hashed.

## 5. Computational closure after the fresh replay

The requested replay is now complete. An independent read-only Julia comparison
of `results/replay_certificate_d7.json` with the preserved original
`results/certificate_d7.json` found **all mathematical fields exactly identical**,
recursively excluding only metadata and timings. This includes all sector error
bounds, partition enclosures, rational defects, and the bound
0.2187941698399417. The recorded replay runtime is 949.34965381 seconds.

The reviewer independently recomputed and checked all 36 source SHA-256 values
in the replay metadata, and both certificate-file hashes recorded in
`results/replay_comparison.json`. All checks passed. The replay log terminates
with the completed six-cell verification and the expected certificate inputs
and bound. The publication workflow also reports 40 passing assertions from
`replay_check.jl`; the independent comparison here was a separate read-only
check and did not rewrite its output or rerun diagonalisation.

The historical source archive is identified by commit
`6c375b5981f1f7c7813bce53ca01a1ce1c04c9f2`; the publication workflow reports this
as an identity-only commit rewrite preserving the historical tree. The fresh
replay uses its own verified source hashes, so its validation does not depend
on that reported identity rewrite.

This closes the pending computational condition stated earlier in this audit.
The delta=7/10 certificate is supported by the mathematical review and fresh
complete replay, subject to the ordinary arithmetic/software assumptions and
non-formal-verification limitations stated above. No mathematical bound drift
or remaining substantive audit blocker was found.
