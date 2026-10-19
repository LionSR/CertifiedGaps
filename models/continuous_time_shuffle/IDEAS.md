# Finite-degree chronological shuffle identities

## Candidate routes

1. Restrict the existing permutation chamber decomposition to permutations preserving each
   bond's internal event order. This proves the identity for arbitrary Banach-valued integrands,
   without imposing commutativity or constructing any physical partition function.
2. Induct on event count using iterated one-dimensional integration. This is mathematically
   direct but makes dependent bounds and Fubini bookkeeping unnecessarily expensive.
3. Use an abstract partition hypothesis. This is too weak: the actual chamber partition and
   within-block order constraints must be proved rather than assumed.

Choose route 1. First prove a general finite ordering-constraint theorem, then exhibit the
measure-preserving coordinate splitting identifying two ordered bond histories with a product
of simplexes. Finite sums grow with the number of shuffles, while Lean verifies the theorem
uniformly in the event counts. No numerical seed or physical trace identity is claimed.

## Outcome

The generic precedence-constraint partition, the two-history product and iterated-integral
identities, and the arbitrary finite-family product-simplex shuffle theorem are proved in
`lean/CertifiedGaps/Transfer/HistoryShuffle.lean`. The finite-family proof constructs a
measure-preserving flattening of event slots and proves the exact within-bond order criterion.
No simplex partition, measure transport, or physical trace identity is assumed.

All event counts, including zero counts and an empty family, are covered uniformly. The
infinite summation is proved separately in `Transfer/ProductPoissonHistory.lean`. Its route
uses the measurable partition by complete word vectors to derive the actual finite-product
measure decomposition, then proves absolute convergence before ungrouping the chronological
series. This avoids assuming that products and countable sums of measures commute. A canonical
zero extension connects each subtype component with the ambient finite-family shuffle theorem.
The resulting series has an explicit exponential majorant. No physical partition-function
expansion is claimed. Compile, axiom-policy, and kernel-replay checks are documented in the
accompanying README.
