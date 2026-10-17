# Method: core pipeline, checker semantics, and the regeneration audit

## Pipeline

Each grid point is one call of the shared core:

```julia
model = RingModel(2, -kron(Z, Z), -g * X - h * Z)          # exact, physical σᶻ frame
cert  = certify_gap(model; n0 = 6, β0 = BETA0[(g, h)], rp = rp_certificate(model; factors = [Z]))
```

`certify_gap` performs these steps, all in `julia/src`:

1. **RP.** `verify_rp` checks $-h_{\text{bond}} = Z \otimes Z$ and the reality of $o = -gX - hZ$
   exactly.
2. **Sectors.** `ring_blocks` builds the translation orbits of the $2^L$ configurations and, for
   each momentum, the block $`\langle s,k|H|r,k\rangle = \sqrt{p_r/p_s}\sum_l c(r \to T^l s)\,e^{ikl}`$
   in interval arithmetic (exact couplings, enclosed phases and norm ratios). Its `Float64` center
   is exactly Hermitian (real for $k = 0, L/2$), and the Frobenius norm of the entrywise radii
   bounds $\|A - M\|_2$. The model is real, so $`A_{-k} = \overline{A_k}`$ has the same spectrum and
   only $0 \le k \le L/2$ is built, with multiplicity 2 for $0 < k < L/2$. No charge is conserved
   for $h > 0$, so there is one sector. At $L = 12$: 7 blocks, the largest of dimension 352.
3. **Enclosure.** `verified_eigen` turns the LAPACK eigenpairs of each center into a guarantee:
   residual $R = MX - X\Lambda$ and orthogonality defect $G = X^\dagger X - I$ in `Float64` with
   the a priori bound $\gamma_K \sum |a_i b_i|$ ($K = N + 8$ real, $K = 2N + 8$ complex through four
   real GEMMs) and Frobenius norms; then the smaller of the Ostrowski and polar-factor radii,
   plus the construction error, by Weyl. Radii at $L = 12$ are about $10^{-9}$.
4. **Seeds.** `log_partition_enclosure` encloses $\log Z_6(\beta_0)$, $\log Z_{12}(\beta_0)$ and
   $\log Z_{12}(2\beta_0)$; the defect enclosures' upper ends are rounded up to the grid
   $10^{-15}$ to give the exact $\varepsilon_{\mathrm{s}}, \varepsilon_\tau$.
5. **Chain.** `rational_gap_chain` (the number of doublings that gives the best bound) yields the
   exact strict lower bound `gap_lower_bound`; `certified_gap_bound` gives the interval bound
   `interval_bound_decimal` ($\Delta_L \ge$ it), used as the admissibility screen. The README
   headline quotes the interval bound and lists the chain bound beside it.

The certificate stores the model, RP data, seed, $\log Z$ intervals, spectrum summaries (block
count, largest block, largest radius) and chain, plus the fields $(g, h)$ and the run metadata
with `model_provenance` hashes.

## Checker semantics

`check.jl` (fast) enforces the model-local source hashes, compares the stored model, RP factor
and seed exactly with `longitudinal(g, h)`, `C = Z` and `BETA0`, and reruns the core
`check_certificate(cert; spectra = false)`: the exact RP identity and the exact chain, which must
reproduce the stored bound from the stored seed. `check.jl --full` adds `spectra = true`: both
seed spectra are re-enclosed from scratch and the fresh defect enclosures must not exceed the
stored $\varepsilon_{\mathrm{s}}, \varepsilon_\tau$. This is the sound replay criterion: the stored
defects are the only spectral input of the chain, so re-proving that they are upper bounds
re-proves the certificate, independent of the LAPACK build that produced it. (The earlier local
checker instead required each fresh eigenvalue interval to lie inside a stored interval of
doubled radius $2\eta$; the core certificate stores no eigenvalues.)

## Seed selection

`BETA0` was selected by the earlier floating-point scan over $`n_0 \in \{2, 4, 6\}`$ and
$`\beta_0 \in \{1/4, \dots, 64\}`$ (best `certified_gap_bound`), kept unchanged so that the seeds of
the Lean statement stay valid. `scripts/scan.jl` repeats the scan from the core momentum blocks
and matches 29 of 30 points; at $g = 1/4, h = 1$ its `Float64` defects are at rounding level and
its preferred $\beta_0 = 23/4$ carries a meaningless float bound. With the tighter enclosures some
points would now admit a better $\beta_0$; that was not pursued.

## Regeneration with the core (old versus new)

Old: the local pipeline (radius $2\eta$ with $\gamma_{16d+32}$ allowances), interval bound
`certified_gap_bound`. New: the core pipeline at the same $(n_0, \beta_0)$; "chain" is the exact
`gap_lower_bound`, "interval" the core `interval_bound_decimal`. Lean $s_0, t_0$ are the seed
defects assumed in `IsingLongitudinal.lean`; every new $\varepsilon \le$ Lean value exactly.

| $g$ | $h$ | old $\varepsilon_{\mathrm{s}}$ | new $\varepsilon_{\mathrm{s}}$ | Lean $s_0$ | old $\varepsilon_\tau$ | new $\varepsilon_\tau$ | Lean $t_0$ | old bound | chain | interval |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 1//4 | 1//20 | 9.064e-06 | 1.214e-06 | 9.064e-06 | 3.14e-05 | 4.17e-08 | 3.141e-05 | 0.3923 | 0.5319 | 0.5324 |
| 1//4 | 1//10 | 4.674e-06 | 1.79e-07 | 4.675e-06 | 1.799e-05 | 2.473e-08 | 1.799e-05 | 0.7332 | 1.0734 | 1.0743 |
| 1//4 | 1//4 | 2.225e-06 | 1.173e-08 | 2.226e-06 | 8.857e-06 | 1.319e-08 | 8.857e-06 | 1.6606 | 2.6111 | 2.6130 |
| 1//4 | 1//2 | 2.063e-06 | 3.287e-09 | 2.063e-06 | 8.242e-06 | 1.31e-08 | 8.242e-06 | 1.9757 | 3.1451 | 3.1474 |
| 1//4 | 1//1 | 2.494e-06 | 4.133e-09 | 2.494e-06 | 9.964e-06 | 1.647e-08 | 9.964e-06 | 2.0336 | 3.2512 | 3.2536 |
| 1//2 | 1//20 | 0.001314 | 0.001309 | 0.001314 | 1.872e-05 | 6.067e-07 | 1.872e-05 | 0.4564 | 0.4567 | 0.4572 |
| 1//2 | 1//10 | 0.0002267 | 0.0002239 | 0.0002267 | 1.128e-05 | 3.349e-08 | 1.129e-05 | 0.9592 | 0.9631 | 0.9639 |
| 1//2 | 1//4 | 7.236e-06 | 5.595e-06 | 7.237e-06 | 7.137e-06 | 5.785e-07 | 7.138e-06 | 2.3616 | 2.5189 | 2.5217 |
| 1//2 | 1//2 | 1.625e-06 | 3.667e-08 | 1.625e-06 | 6.516e-06 | 1.707e-07 | 6.517e-06 | 2.7754 | 3.6860 | 3.6891 |
| 1//2 | 1//1 | 1.974e-06 | 3.535e-09 | 1.974e-06 | 7.885e-06 | 1.418e-08 | 7.885e-06 | 2.7276 | 4.3048 | 4.3080 |
| 3//4 | 1//20 | 0.0371 | 0.03709 | 0.0371 | 0.0003653 | 0.000353 | 0.0003654 | 0.3170 | 0.3167 | 0.3171 |
| 3//4 | 1//10 | 0.008573 | 0.00857 | 0.008573 | 6.931e-05 | 6.112e-05 | 6.931e-05 | 0.7712 | 0.7708 | 0.7714 |
| 3//4 | 1//4 | 0.0002371 | 0.0002357 | 0.0002372 | 0.000133 | 0.0001274 | 0.0001331 | 2.1018 | 2.1048 | 2.1066 |
| 3//4 | 1//2 | 4.347e-06 | 2.835e-06 | 4.348e-06 | 1.065e-05 | 4.605e-06 | 1.065e-05 | 3.0117 | 3.2174 | 3.2211 |
| 3//4 | 1//1 | 1.856e-06 | 3.028e-08 | 1.856e-06 | 7.341e-06 | 4.862e-08 | 7.342e-06 | 3.1374 | 4.5166 | 4.5202 |
| 1//1 | 1//20 | 0.184 | 0.184 | 0.1841 | 0.007424 | 0.007413 | 0.007425 | 0.0763 | 0.0763 | 0.0763 |
| 1//1 | 1//10 | 0.05123 | 0.05122 | 0.05123 | 0.004667 | 0.004659 | 0.004667 | 0.5043 | 0.5039 | 0.5043 |
| 1//1 | 1//4 | 0.002469 | 0.002467 | 0.002469 | 0.002244 | 0.002239 | 0.002245 | 1.6092 | 1.6076 | 1.6097 |
| 1//1 | 1//2 | 5.079e-05 | 4.936e-05 | 5.079e-05 | 9.929e-05 | 9.358e-05 | 9.929e-05 | 2.7523 | 2.7672 | 2.7693 |
| 1//1 | 1//1 | 2.206e-06 | 5.276e-07 | 2.207e-06 | 7.689e-06 | 9.832e-07 | 7.69e-06 | 3.6403 | 4.2822 | 4.2864 |
| 3//2 | 1//20 | 0.1477 | 0.1477 | 0.1478 | 0.05423 | 0.05423 | 0.05424 | 0.1692 | 0.1690 | 0.1692 |
| 3//2 | 1//10 | 0.09847 | 0.09847 | 0.09848 | 0.04084 | 0.04083 | 0.04085 | 0.3601 | 0.3596 | 0.3601 |
| 3//2 | 1//4 | 0.01818 | 0.01818 | 0.01818 | 0.01277 | 0.01277 | 0.01278 | 1.0556 | 1.0546 | 1.0557 |
| 3//2 | 1//2 | 0.001195 | 0.001193 | 0.001195 | 0.001672 | 0.001666 | 0.001673 | 2.1293 | 2.1278 | 2.1303 |
| 3//2 | 1//1 | 2.24e-05 | 2.076e-05 | 2.241e-05 | 3.344e-05 | 2.689e-05 | 3.345e-05 | 3.7050 | 3.7655 | 3.7706 |
| 2//1 | 1//20 | 0.03251 | 0.03251 | 0.03252 | 0.03673 | 0.03672 | 0.03673 | 0.9121 | 0.9114 | 0.9121 |
| 2//1 | 1//10 | 0.0294 | 0.0294 | 0.02941 | 0.03285 | 0.03284 | 0.03285 | 0.9643 | 0.9636 | 0.9643 |
| 2//1 | 1//4 | 0.01577 | 0.01577 | 0.01578 | 0.01668 | 0.01668 | 0.01669 | 1.2728 | 1.2715 | 1.2729 |
| 2//1 | 1//2 | 0.003329 | 0.003327 | 0.003329 | 0.003195 | 0.003189 | 0.003196 | 2.0005 | 1.9996 | 2.0010 |
| 2//1 | 1//1 | 0.0001302 | 0.0001286 | 0.0001302 | 0.0003861 | 0.0003799 | 0.0003862 | 3.4912 | 3.4957 | 3.4989 |

All 60 defects decreased strictly and all 30 interval bounds increased. The exact chain bound is
below the old interval bound at 13 points, by at most $`0.13\%`$ (relative); it is a separate,
slightly more conservative exact rational evaluation of the recurrence, and it is the bound the
certificate claims.
