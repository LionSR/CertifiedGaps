# Candidate routes (written before implementation)

Task 4 tests the proposed Lean seed certificate; it does not claim a completed gap proof.
All new files remain in this directory; the original Lean project is read-only.

1. Full base-3 packing: digit-pattern masks have closed geometric-series recipes and
   nearest-neighbor transitions become fixed shifts. Cost is O(n) bignum operations per
   matvec on 3^n lanes. Chosen as the first executable prototype.
2. Compact magnetization packing: reduces n=12 from 531441 to 73789 lanes, but rank
   compression destroys constant shift distances. The advertised O(50) shifts needs an
   additional construction; audit this rather than assuming it.
3. Nonnegative shifted/Marshall-transformed operator: removes signed intermediate
   subtraction and permits cheap power-filter generation from product states. Signed
   arithmetic still needs explicit lane bounds; use separate positive/negative limbs.
4. Ritz plus residual: prove exact integer Gram inequalities first. Global low-eigenvalue
   identification additionally requires a counted anchor or a complementary-space bound.
   A residual alone is not a lower bound on the lowest eigenvalues.
5. Full spectral certification: exact LDL/inertia at rational endpoints is a fallback,
   but sparse elimination fill and kernel cost may defeat the proposed schedule.
6. Minimal honest baseline: kernel-check H6/H8 arithmetic and crude global bounds,
   recording explicitly any missing bridge or strong low-level lower certificate.

Scaling variables: vector length 3^n, lane width set by filter growth and dot-product
carry bounds, k retained vectors, filter degree, and recursive anchor count. Beta affects
required spectral precision and energy cutoff, not a single matvec directly.

Experiments using Float64 are non-certifying. Julia is only a generator; no generated
number is accepted as a theorem without a Lean kernel proof. No native_decide or sorry.

## Outcome

Route 1 (full base-3 packing) works as an arithmetic primitive: dense n = 12 matvec + dot
checks by `decide +kernel` in tens of seconds, and power-filter Ritz/residual witnesses for
H6, H8, H10 are kernel-checked. Route 2 fails as stated: compact Sz = 0 packing needs 3804
distinct shifts at n = 12, not about 50. Route 4 is blocked by the anchor: with the untuned
cut comparison (two open half chains minus 4I), a Float64 scan of the full H6 and H8 spectra
finds no rank k with the separation the counting lemma needs, even for exact Ritz vectors.
Routes 5 and 6 were not completed. No seed, eigenvalue enclosure or gap is certified; the
long review and the implementation budget (8.5–15k Lean lines) are in METHOD.md.
