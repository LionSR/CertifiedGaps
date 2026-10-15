# Method: the two-level route, the all-model inventory and the remaining proof work

This is the long companion to [README.md](README.md). Everything below that is not a Lean
theorem is labelled; the floating-point numbers are planning evidence only.

## A cheaper spectral route for the selected point

The selected point is the longitudinal Ising ring `H = −Σ ZZ − g Σ X − h Σ Z` at
`g = 1/4, h = 1/20`, seed `(n₀, β₀) = (6, 97/4)`. Two approximate eigenvectors per ring plus a
classical diagonal comparison suffice numerically; full Rayleigh–Ritz matrices, a recursive
half-chain anchor, momentum sectors and SU(2) are unnecessary at this point.

**Classical anchor.** Let `D` be the classical diagonal part of `H`. Since `−gX ≥ −gI` on one
spin, `H ≥ D − ngI`. Exactly two periodic configurations (all up, all down) have no domain
wall; every other one has at least two, so its classical energy is at least
`−n(1+h) + 4 + 2h`. At `g = 1/4, h = 1/20`:

```
ρₙ = −n(1+h+g) + 4 + 2h = (82 − 26n)/20,     λ₃(H₆) ≥ −3.7,     λ₃(H₁₂) ≥ −11.5.
```

The direct classical counts (`anchor6`, `anchor12`) are proved in Lean. The operator
inequality and the eigenvalue-count consequence (min–max) are still proof obligations.

**Residual intervals.** For the integer Hamiltonian `A = 20H` and scale `S = 2⁵²`, the kernel
checks

```
|vᵢᵀvᵢ − S²| ≤ 10⁻¹² S²,   |v₁ᵀv₂| ≤ 10⁻¹² S²,   Σᵢ ‖S A vᵢ − λᵢ vᵢ‖² ≤ 10⁻¹⁸ S⁴,
λ₁ + S < λ₂,               λ₂ + S < (20ρ) S.
```

Once the norm identities are proved, each residual gives an eigenvalue within
`10⁻⁹/(20√(1−10⁻¹²)) < 5.01×10⁻¹¹` of the candidate `μᵢ = λᵢ/(20S)`. The two intervals are
disjoint and lie below ρ, so a proved count of at most two eigenvalues below ρ makes them
exactly the first two. Cross-Gram orthogonality is then unnecessary; it is checked only for
comparison. For larger spectra the analogous saving is to certify small clusters of
overlapping intervals and their multiplicities instead of a full k×k Gram matrix; that still
needs a counting theorem and a reliable anchor.

**Trace sandwich.** With a deliberately conservative energy radius `r = 10⁻¹⁰`,

```
Σᵢ₌₁² exp(−β(μᵢ + r)) ≤ Zₙ(β) ≤ Σᵢ₌₁² exp(−β(μᵢ − r)) + (2ⁿ − 2) exp(−βρₙ).
```

`scripts/feasibility.jl` evaluates this in Float64 (**non-certifying**; it assumes the spectral
theorem above), writing `results/feasibility.json`:

| Quantity | Sandwich result | Required target |
|---|---:|---:|
| Spatial defect upper bound | 1.210665×10⁻⁶ | 9.064×10⁻⁶ |
| Thermal defect upper bound | 9.949076×10⁻⁹ | 3.141×10⁻⁵ |
| 6-site normalised tail at β | 2.7771×10⁻²⁷ | included above |
| 12-site normalised tail at β | 1.2424×10⁻¹⁰ | included above |

The predicted log-partition intervals also fit strictly inside the four endpoints that
`LogEnclosures` asks for, so what remains is a proof task, not a need for sharper vectors.
The `h = 1/10` point has only one level below the inventory cutoff, but the rank-two classical
anchor does not isolate that level by itself; `h = 1/20` is the cheapest point found with a
complete numerical anchor argument.

**Scalar exponentials.** A concrete route is `Real.exp_bound` on each exponent divided by 1024
(the largest exponent here is about 620), followed by ten outward-rounded rational squarings.
A 16-term Taylor proposal is amply sufficient at these margins, but its rounding and remainder
checks must themselves be proved.

## Inventory of all 43 seed targets

`scripts/inventory.jl` reads the stored spectra of every certified model (30 longitudinal
Ising points, 5 staggered-Heisenberg, 3 spin-½ dimerised, 2 large-D, and one each of the
Heisenberg ladder, spin-1 dimerised δ = 9/10 and Ising ladder g = 3) and writes
`results/inventory.{json,tsv}`. The seed defects are not eigenvalue tolerances; to avoid
silently using all of a seed's slack it sets

```
σ = min(εₛ − εₛ,source, εₜ − εₜ,source),   η = min(10⁻⁶, σ/16),   r = η/(4β),
```

and counts the levels below `E₀ + (log D + log(1/η))/β`. "source" is the stored external
enclosure, not a Lean theorem. Counts use the nominal stored spectra with multiplicities and
are **non-certifying**.

| Family | Largest full dimension | Largest Sz block, when used | k at the cutoff | Full three-Gram arithmetic scenario |
|---|---:|---:|---:|---:|
| Longitudinal Ising, 30 points | 4,096 | 4,096 (no Sz reduction) | 1–145 | 0.014–150 s per point |
| Spin-1 dimerised δ = 0.9 | 6,561 | 1,107 | 195 | 5–80 s |
| Staggered Heisenberg, 5 points | 65,536 | 12,870 | 424–1,249 | 0.14–19.8 h |
| Heisenberg ladder J = 4 | 65,536 | 12,870 | 1,238 | 1.16–19.4 h |
| Spin-½ dimerised, 3 points | 65,536 | 12,870 | 862–1,558 | 0.56–30.8 h |
| Ising ladder g = 3 | 65,536 | 65,536 (Sz not conserved) | 322 | 15.7 h |
| Large D = 4 | 531,441 | 73,789 | 1,480 | 15.4–379 h |
| Large D = 3.206 | 531,441 | 73,789 | 2,465 | 42.5–1,051 h |

These are **arithmetic-only scenario estimates** (`scripts/costs.jl`), not measured certificate
run times or proven lower bounds. For conserved charge they range from k spread evenly over
the charge sectors to all k in the largest block, pricing each product at the largest block
width. They exclude vector generation, orthogonalisation, packing, anchors, spectral proofs
and scalar exponential bounds; no momentum, cyclotomic or SU(2) saving is assumed. The low β
of the large-D seeds is why thousands of levels matter despite small seed lengths. These
numbers reject a blanket 15–60 minute estimate for these models with full k×k Gram checks;
they do not refute a separate spin-1 Heisenberg proposal at β = 49/4.

### Kernel multiplication benchmark

`scripts/costs.jl` checks `(2ᵇ−1)(2ᵇ−3) = 2²ᵇ − 2ᵇ⁺² + 3` by `decide +kernel` on dense
operands (Linux run recorded in `results/costs.json`):

| Operand bits | Kernel time |
|---:|---:|
| 262,144 | 0.00111 s |
| 1,179,648 | 0.00471 s |
| 3,706,560 | 0.0304 s |
| 21,251,232 | 0.415 s |

These are primitive microbenchmarks; `costs.jl` extrapolates with a local size exponent 1.15
from the nearest measured size.

### Full table

`julia --project=. scripts/report_tables.jl` prints this table from `results/costs.json`.
`η` and `r` as above; Gram minutes exclude vector generation, anchors, analytic proofs and
encoding costs.

| Model / parameter | D | β | εₛ | εₜ | η | r | cutoff − E₀ | k | Gram minutes, balanced–concentrated |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| LargeD 4_1 | 531441 | 2.25 | 0.06411 | 0.04053 | 3.79e-08 | 4.21e-09 | 13.454 | 1480 | 924–2.27e+04 |
| LargeD 1603_500 | 531441 | 3.048 | 0.19164415 | 0.046462453 | 3.66e-11 | 3e-12 | 12.21 | 2465 | 2.55e+03–6.31e+04 |
| Heisenberg ladder J4 | 65536 | 2 | 0.03509 | 0.02463 | 3.63e-07 | 4.53e-08 | 12.96 | 1238 | 69.5–1.17e+03 |
| Staggered h1_1 | 65536 | 5.25 | 0.0006937 | 0.0007327 | 8.26e-10 | 3.94e-11 | 6.096 | 424 | 8.36–137 |
| Staggered h3_4 | 65536 | 5.25 | 0.002343 | 0.003123 | 3e-08 | 1.43e-09 | 5.4117 | 534 | 13.2–217 |
| Staggered h1_2 | 65536 | 5.5 | 0.01075 | 0.01116 | 1.87e-07 | 8.51e-09 | 4.8331 | 786 | 28.2–470 |
| Staggered h7_20 | 65536 | 5.75 | 0.03266 | 0.02744 | 2.56e-07 | 1.11e-08 | 4.5686 | 1152 | 60.2–1.01e+03 |
| Staggered h1_4 | 65536 | 6.25 | 0.07987 | 0.04477 | 2.22e-07 | 8.87e-09 | 4.2259 | 1249 | 70.7–1.19e+03 |
| Dimerized half d1_2 | 65536 | 6 | 0.01133 | 0.01224 | 1.94e-07 | 8.07e-09 | 4.4247 | 862 | 33.9–565 |
| Dimerized half d3_10 | 65536 | 6.5 | 0.09516 | 0.0518 | 5.3e-08 | 2.04e-09 | 4.2837 | 1558 | 110–1.85e+03 |
| Dimerized half d1_4 | 65536 | 7.25 | 0.1658 | 0.05123 | 1.32e-07 | 4.55e-09 | 3.7146 | 913 | 38–634 |
| Dimerized one d9/10 | 6561 | 3.5 | 0.04972 | 0.03206 | 2.71e-07 | 1.94e-08 | 6.8315 | 195 | 0.0845–1.33 |
| Ising ladder G3 | 65536 | 1.5 | 0.1125 | 0.05162 | 2.98e-07 | 4.96e-08 | 17.412 | 322 | 941–941 |
| Ising g1/4 h1/20 | 4096 | 24.25 | 9.064e-06 | 3.141e-05 | 1.41e-11 | 1.45e-13 | 1.3734 | 2 | 0.000706–0.000706 |
| Ising g1/4 h1/10 | 4096 | 13.75 | 4.675e-06 | 1.799e-05 | 3.47e-11 | 6.31e-13 | 2.3565 | 1 | 0.000235–0.000235 |
| Ising g1/4 h1/4 | 4096 | 6.5 | 2.226e-06 | 8.857e-06 | 1.95e-11 | 7.5e-13 | 5.0736 | 25 | 0.0765–0.0765 |
| Ising g1/4 h1/2 | 4096 | 5.5 | 2.063e-06 | 8.242e-06 | 4.56e-12 | 2.07e-13 | 6.2602 | 25 | 0.0765–0.0765 |
| Ising g1/4 h1/1 | 4096 | 5.25 | 2.494e-06 | 9.964e-06 | 3.56e-12 | 1.7e-13 | 6.6053 | 13 | 0.0214–0.0214 |
| Ising g1/2 h1/20 | 4096 | 13 | 0.001314 | 1.872e-05 | 3.82e-12 | 7.35e-14 | 2.6622 | 3 | 0.00141–0.00141 |
| Ising g1/2 h1/10 | 4096 | 8 | 0.0002267 | 1.129e-05 | 5.98e-10 | 1.87e-11 | 3.6945 | 12 | 0.0184–0.0184 |
| Ising g1/2 h1/4 | 4096 | 4.5 | 7.237e-06 | 7.138e-06 | 3.25e-11 | 1.8e-12 | 7.2153 | 74 | 0.654–0.654 |
| Ising g1/2 h1/2 | 4096 | 4 | 1.625e-06 | 6.517e-06 | 1.48e-12 | 9.24e-14 | 8.8894 | 49 | 0.288–0.288 |
| Ising g1/2 h1/1 | 4096 | 4 | 1.974e-06 | 7.885e-06 | 7.46e-12 | 4.66e-13 | 8.4847 | 25 | 0.0765–0.0765 |
| Ising g3/4 h1/20 | 4096 | 8 | 0.0371 | 0.0003654 | 4.67e-09 | 1.46e-10 | 3.4374 | 21 | 0.0544–0.0544 |
| Ising g3/4 h1/10 | 4096 | 5.25 | 0.008573 | 6.931e-05 | 2.96e-10 | 1.41e-11 | 5.7632 | 72 | 0.619–0.619 |
| Ising g3/4 h1/4 | 4096 | 3.5 | 0.0002372 | 0.0001331 | 3.16e-09 | 2.26e-10 | 7.9685 | 103 | 1.26–1.26 |
| Ising g3/4 h1/2 | 4096 | 3.5 | 4.348e-06 | 1.065e-05 | 4.52e-11 | 3.23e-12 | 9.1823 | 68 | 0.552–0.552 |
| Ising g3/4 h1/1 | 4096 | 3.5 | 1.856e-06 | 7.342e-06 | 2.56e-11 | 1.83e-12 | 9.3445 | 25 | 0.0765–0.0765 |
| Ising g1/1 h1/20 | 4096 | 6.5 | 0.1841 | 0.007425 | 5.97e-08 | 2.3e-09 | 3.8387 | 25 | 0.0765–0.0765 |
| Ising g1/1 h1/10 | 4096 | 4.25 | 0.05123 | 0.004667 | 2.98e-08 | 1.75e-09 | 6.0344 | 75 | 0.671–0.671 |
| Ising g1/1 h1/4 | 4096 | 3 | 0.002469 | 0.002245 | 2.01e-08 | 1.68e-09 | 8.68 | 145 | 2.49–2.49 |
| Ising g1/1 h1/2 | 4096 | 3 | 5.079e-05 | 9.929e-05 | 1.1e-10 | 9.17e-12 | 10.416 | 122 | 1.77–1.77 |
| Ising g1/1 h1/1 | 4096 | 3 | 2.207e-06 | 7.69e-06 | 5.42e-11 | 4.51e-12 | 10.652 | 37 | 0.166–0.166 |
| Ising g3/2 h1/20 | 4096 | 3.5 | 0.1478 | 0.05424 | 4.36e-07 | 3.11e-08 | 6.5612 | 56 | 0.376–0.376 |
| Ising g3/2 h1/10 | 4096 | 3.25 | 0.09848 | 0.04085 | 5.67e-07 | 4.36e-08 | 6.985 | 60 | 0.431–0.431 |
| Ising g3/2 h1/4 | 4096 | 2.75 | 0.01818 | 0.01278 | 1.55e-07 | 1.41e-08 | 8.7271 | 92 | 1.01–1.01 |
| Ising g3/2 h1/2 | 4096 | 2.5 | 0.001195 | 0.001673 | 1.52e-08 | 1.52e-09 | 10.529 | 93 | 1.03–1.03 |
| Ising g3/2 h1/1 | 4096 | 2.5 | 2.241e-05 | 3.345e-05 | 4.57e-10 | 4.57e-11 | 11.93 | 61 | 0.445–0.445 |
| Ising g2/1 h1/20 | 4096 | 2.25 | 0.03252 | 0.03673 | 1.73e-07 | 1.92e-08 | 10.618 | 110 | 1.44–1.44 |
| Ising g2/1 h1/10 | 4096 | 2.25 | 0.02941 | 0.03285 | 2.93e-07 | 3.25e-08 | 10.383 | 92 | 1.01–1.01 |
| Ising g2/1 h1/4 | 4096 | 2.25 | 0.01578 | 0.01669 | 3.22e-07 | 3.58e-08 | 10.341 | 81 | 0.782–0.782 |
| Ising g2/1 h1/2 | 4096 | 2.25 | 0.003329 | 0.003196 | 1.57e-08 | 1.75e-09 | 11.683 | 85 | 0.861–0.861 |
| Ising g2/1 h1/1 | 4096 | 2 | 0.0001302 | 0.0003862 | 2.91e-09 | 3.64e-10 | 13.986 | 79 | 0.744–0.744 |

## Implementation observations

- A single-expression direct 12-site checker was stopped at 7.4 GB RSS and 38 CPU s without
  completing. Splitting the scalar work into 64 independently checked blocks
  (`GeneratedChunks`) succeeds; `ChunkCoverage` proves the block list is exactly the fixed
  range enumeration, so a generator cannot duplicate or omit a block.
- The first packed version checked the vector arithmetic in about 0.4 s but spent about 140 s
  elaborating two padded decimal literals. Storing compact 56-bit slots and widening them to
  288-bit lanes inside Lean (`Packed.widen`) cut elaboration to about 6 s.
- Do not treat packed primitive timings as full-checker timings: vector preparation, Gram
  products, classical counts, scalar bounds and peak kernel memory are separate costs.

## Remaining proof work for the selected point

The missing theorem is `LogEnclosures` (README, "Seeds and certified bounds"). A route to it:

| Component | Additional Lean lines (estimate) | Obligation |
|---|---:|---|
| Matrix identification and scalar-sum interpretation | 200–400 | Binary states ↔ `Fin n → Fin 2`; prove `A = 20H` against `RPRing.hamiltonian`; identify the checked sums as norms and Gram entries |
| Residual-to-spectrum bridge | 150–350 | Apply `exists_near_of_weighted_residual` through the Hermitian spectral theorem |
| Classical comparison and counting | 250–500 | Prove `H ≥ D − ngI`, transfer the classical count by min–max, exhaust two eigenvalues |
| Rational exponential and trace bounds | 200–450 | Outward scalar exp bounds and the trace sandwich; derive the four log bounds |
| Integration | 100–250 | Connect to `IsingLongitudinal.gap` without the seed `sorry` |
| Optional packed-backend correctness | 300–700 | Mask/flip/reverse, carries, middle-product extraction, widening; unnecessary with the chunk route |

That is an estimated 900–1,950 lines by the direct route (1,200–2,650 with a reusable packed
backend), not a measured completion. A backend covering all seven model families needs
symmetry identification, anchors, cluster multiplicities and scalar bounds; 8–15k lines is the
more defensible planning figure. The larger Gram scenarios rule out a uniform minute-scale
run time without further structural reductions.

## Suggested order

1. The selected longitudinal Ising point: explicit two-state anchor and a large seed margin.
   A separate theorem for this point can become unconditional before the 30-point `seeds`.
2. The other longitudinal Ising points and spin-1 dimerisation δ = 0.9 (largest Sz block 1,107).
   General points still need anchors; the spin-1 Sz basis change is extra work.
3. Strong staggered fields (h = 1, 3/4), then the Heisenberg ladder and the easier spin-½
   dimerisation targets, prioritised by cluster counts and anchor strength.
4. Weaker staggered fields, harder dimerisation targets and the Ising ladder (no Sz
   decomposition; the most expensive Gram scenario among the spin-½ models).
5. Large D = 4, then D = 3.206; the latter has only about 5.85×10⁻¹⁰ slack against its stored
   enclosure, so redesigning its seed is more promising than a dense Gram plan.
