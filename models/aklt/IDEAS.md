# Routes considered

This memo was reconstructed after the fact from the README, the 2+1D note
(`notes/thermal-mpo-gap-2plus1d.tex`) and the history of this directory. The pilot predates the
rule that every model directory starts with a route memo.

The target is the spin-3/2 AKLT model on the honeycomb lattice in projector normalization,
\(H=\sum_{\langle ij\rangle}P_{ij}^{(3)}\), with exact ground energy 0. The 2+1D note replaces the
two purities by the partition-function excess \(R=Z-1\). Temporal doubling squares \(R\) from
above, and two scalar spatial gluing inequalities would give \(R'\le5R^2\) for \(R\le0.01\). A
uniform gap then needs one certified finite seed and a proof of the gluing inequalities.

* **Reflection positivity of the bond (the 1D criterion's route to gluing).** Fails for the
  standard construction. With time reversal \(\theta(A)=U\bar AU^\dagger\),
  \(U=e^{-i\pi S^y}\), and the quadrupole \(Q=(S^z)^2-\tfrac54\), the note shows
  \(\operatorname{Tr}[e^{-tP^{(3)}}(Q\otimes\theta Q)]=\tfrac75(e^{-t}-1)<0\). This rules out that
  local RP construction, not the scalar gluing inequalities themselves.
* **Self-adjoint spatial transfer for a projection filter.** Checked at width two (a transfer of
  dimension 256) with local Python helpers. No positive metric makes \(T\) or \(T^2\)
  self-adjoint, so the self-adjoint route fails for that filter. Exact characteristic polynomials
  still give \(d_{2L}\le d_L^2/4\) at width 2 for every even \(L\ge4\), but this does not extend to
  the mixed-width, mixed-depth tuples the recursion needs.
* **Exact diagonalization for a seed.** The 8-spin (2×2-cell) quotient has 65,536 states and is
  cheap by total-\(S^z\) sectors. 16 spins (\(4^{16}\) states) is at the edge of dense ED and
  32 spins is out of reach. Small periodic quotients have short winding cycles, so they give
  finite-size checks, not bulk evidence.
* **Thermal MPO purification.** The cost is polynomial in the bond dimension χ and linear in
  β/δβ, so it could reach the 16- and 32-spin tori if compression is accurate. Each bond gate is
  exact (\(e^{-tP}=I+\operatorname{expm1}(-t)P\)), so only Trotter and truncation errors remain,
  but this pipeline has no rigorous truncation bound.
* **PEPS boundary methods.** Relative factorization of boundary operators controls ground-space
  projectors, but those estimates are not thermal-purity certificates.

Chosen route: a floating-point pilot. Compute reference full spectra by symmetry sectors on the
8-spin quotient, and run thermal-MPO purification on 8, 16 and 32 spins, to measure the
diagnostics \(R\), \(C\) and \(G\) and to test whether the purification reproduces the reference.

## Outcome

Negative. The purification does not reproduce the 8-spin reference. At β = 2 the exact
\(\log Z=5.31205522615\) is missed by −1.773, −1.424, −0.587 and −0.333 at χ = 16, 32, 64 and 128.
Halving δβ does not help, so the error is compression bias. The gluing diagnostics are therefore
unresolved. The 8-spin spectrum gives a floating-point seed candidate (\(R(64)\approx0.0083\),
gap ≈ 0.0920), but no gluing inequality is proved, so this route cannot produce a certificate.
The next requirement is a contraction method that reproduces the 8-spin reference and keeps that
accuracy as the width grows.
