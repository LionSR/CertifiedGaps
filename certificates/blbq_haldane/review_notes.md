# Independent mathematical and implementation review

The realignment check in `exact_rp` was independently strengthened: writing K for
transpose on the nine-dimensional operator space and v=vec(I), the rational identities
are R(A)=I-K and R(B)=(I+K)/2+v v'. Thus the eigenvalues of R(-h) are 2c-s
(three skew directions), -s (five symmetric-traceless directions), and -4s
(the scalar direction). This proves precisely the asserted coefficient sign conditions.
The quadrupoles must be symmetrized; unsymmetrized S_a S_b would double-count the
antisymmetric contribution when paired with c-s/2.

The highest-weight projector is derived from the exact SU(2) Casimir spectrum, not from
an iterative eigensolver. Its polynomial identity is checked with Rational{BigInt} on
small charge spaces. Numerical projected bases are tested against unreduced complete
spectra, with explicit multiplet weights. Projector roundoff, randomized spanning and
orthonormalization are not rigorous eigenvalue enclosures. The numerical classification
of any resulting partition functions is unchanged.

The earlier extrapolations are not controlled asymptotics. On n=4,6,8, exponential fits
cannot distinguish a massive tail from a preasymptotic algebraic regime. For theta<=-0.5,
the two gap fits even disagree in sign; the required scale is unresolved. A positive
conditional forecast is not evidence of a true seed, and certainly not a certificate.

The theorem requires even seed length in one-site cells. A complete L=14 spectrum adds
an odd-n=7 diagnostic; it does not extend the admissible seed search beyond n=6. The
next eligible seed, n=8, would require complete L=16 spectra in this workflow.

The larger-ring reference routine searches M=0, which contains one state per integer-spin
multiplet. Its recorded `gap` is therefore a multiplet spacing. Identifying it with the
full spectral gap requires a unique singlet ground state; M=0 alone cannot count a
possible SU(2) ground multiplet. This is explicitly stated in the report, and complete
weighted spectra take precedence wherever available. The extrapolations inherit this
uncertainty and do not establish uniqueness or a uniform gap.
