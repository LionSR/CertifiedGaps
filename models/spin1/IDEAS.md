# Routes considered

This memo was reconstructed after the fact from the README, `METHOD.md`, and the history of this
directory. The reproduction predates the rule that every model directory starts with a route memo.

The task was different from a new certificate. **The periodic spin-one Haldane gap** (OpenAI,
September 24, 2026) already proves a uniform gap for the spin-1 Heisenberg ring from two finite
seeds (source 1 at \((n_0,\beta_0)=(2304,784)\), source 2 at \((60,105/4)\)). The goal here was
to recompute every finite input independently, in exact arithmetic, check each integer against
the paper's evidence files, and then see whether the same inputs support stronger constants.

* **Run the upstream Python verifier.** This is cheap and faithful, but it is not independent,
  and the repository has no Python. Rejected as the main route; the upstream files are fetched
  only as reference evidence (`tools/fetch_upstream.jl`).
* **Port the integer Poisson trace recursion to Julia.** The seeds rest on twisted-ring traces at
  lengths 4–12 (source 1) and 6–10 (source 2) with explicit floor-division rounding bounds.
  Exact `Int64`/`Int128`/`BigInt` arithmetic reproduces every integer, so the comparison with the
  evidence is bit-for-bit. The cost is dominated by the source-1 sectors (about 37 minutes on
  4 threads). Chosen. Source 2 uses an unpacked Eisenstein-basis implementation instead of
  the upstream packed lanes, so it is a second, independent implementation.
* **Exact MPS trial contractions.** These give the variational energy shift
  \(a=700741/500000\) at lengths 72 and 120. Charge blocks (largest 162) in `BigInt` cost about
  33 s. An exact SU(2) multiplicity reduction (largest block 56) brings this to 0.09 s and is
  verified by exact tensor relations (see `METHOD.md`). Chosen.
* **Floating-point controls.** ITensors physical MPOs, the integer trial imported as an MPS, QN
  DMRG and a thermal MPO. These are kept as labeled numerical cross-checks, never as certificate
  inputs.
* **Stronger bounds from the same inputs (`experiments/`).** These reuse the frozen finite
  inputs with no new contraction. One route optimizes the scalar polynomial cutoffs and keeps
  unequal defects through the coupled updates. The other lifts the temperature before extracting
  purity bounds, using the length-120 trial and the threefold sector multiplicity. All
  inequalities are decided by rational logarithm and power enclosures.

## Outcome

Positive. Every thermal and trial integer of both sources is reproduced exactly and agrees with
the upstream evidence (27 thermal certificates, two integer MPS trials, both purity seeds and the
bootstrap constants; stored-output audit in `results/verification.json`). The published bounds
\(\Delta_L\ge\log(20)/784\) (even \(L\ge2304\)) and \((4/105)\log(80/79)\) (even \(L\ge60\)) are
therefore reproduced. The experiments strengthen them, from the same finite inputs, to
\(\Delta_L>0.072901771\) for every even \(L\ge36\) (proved in
`experiments/optimized-spin1-bounds.tex` and independently reviewed in #49). The
reflection-positivity input and the analytic transfer and induction arguments are formalized in
`lean/` and derived in `notes/`.
