# Method: closed-form seeds for the transverse-field Ising ring

## The closed form

Write $H_L = -\sum_j Z_j Z_{j+1} - g\sum_j X_j$ on a ring of $L$ sites. With
$X_j = 1 - 2c_j^\dagger c_j$ and the Jordan–Wigner string built from the $X$'s, the bulk bonds
$Z_j Z_{j+1}$ become fermion bilinears. The boundary bond $Z_L Z_1$ picks up the fermion parity
$P = \prod_j X_j$, so $H_L = P_+ H^{\mathrm{NS}} P_+ + P_- H^{\mathrm{R}} P_-$ with
$P_\pm = (1 \pm P)/2$. $H^{\mathrm{NS}}$ is quadratic with antiperiodic fermions
($k = 2\pi(m + 1/2)/L$) and $H^{\mathrm{R}}$ with periodic ones ($k = 2\pi m/L$). A Bogoliubov
rotation diagonalizes each as $\sum_k \varepsilon_k(\eta_k^\dagger\eta_k - 1/2)$ with
$`\varepsilon_k = 2\sqrt{1 + g^2 - 2g\cos k}`$, except for the unpaired Ramond modes $k = 0$ and
$k = \pi$, which are already diagonal with energies $2(g - 1)$ and $2(g + 1)$ (signed). The
rotation preserves the fermion parity, and with these signed energies the Bogoliubov vacuum is
even in each sector (the ED check below fixes these conventions). Hence

```math
\mathrm{Tr}\, P_\pm e^{-\beta H^{\mathrm{A}}} = \tfrac12\Big[\prod_k 2\cosh\tfrac{\beta\varepsilon_k}{2} \pm \prod_k 2\sinh\tfrac{\beta\varepsilon_k}{2}\Big],
```

and $Z_L = \mathrm{Tr}\,P_+ e^{-\beta H^{\mathrm{NS}}} + \mathrm{Tr}\,P_- e^{-\beta H^{\mathrm{R}}}$ is the formula
in `README.md` (Lieb–Schultz–Mattis 1961, Pfeuty 1970). For $g > 1$ the signed zero-mode energy
$2(g - 1)$ is positive, so all $\varepsilon_k$ are given by the square root. The signs and sectors
were checked against dense exact diagonalization. The tests do this for $L = 4, \dots, 10$ at
$g \in \{11/10, 5/4, 2, 3\}$ and three temperatures, and also check that the ground energy is
$-\frac12\sum_{\mathrm{NS}}\varepsilon_k$. With the wrong zero-mode sign the formula already fails at
$L = 4$ by $`0.65\text{–}68\%`$ (relative error in $Z_4$; prototype run, not part of the tests).

## Evaluation

`logZ` encloses each $\varepsilon_k$ (MPFR $\pi$, $\cos$, $\sqrt{}$ with directed rounding), then
$\log 2\cosh(x/2) = x/2 + \log(1 + e^{-x})$ and $\log 2\sinh(x/2) = x/2 + \log(-\mathrm{expm1}(-x))$.
It sums these into $C, S$ (NS) and $C_R, S_R$ (R) and returns

```math
\log Z = C + \log\tfrac12\big(1 + e^{S - C} + e^{C_R - C}\,(-\mathrm{expm1}(S_R - C_R))\big).
```

All terms are positive for $g > 1$, so there is no cancellation. At 512 bits the enclosures have
width about $10^{-150}$, while the seed defects are $10^{-13}$ to $10^{-38}$.
`Seeds.seed_defects` (core) turns the three enclosures into defect enclosures. Their upper ends
are rounded up to 8 significant digits.

## Two BigFloat precision traps

Both were found while writing the Arb cross-check, and both matter only for `BigFloat` data
whose precision differs from the global one:

1. `Rational{BigInt}(x::BigFloat)` rounds to the global precision, not to `precision(x)`. For
   example, a 512-bit $1/3$ converted at the default 256 bits differs from its exact value by
   $1.4 \times 10^{-78}$. `exact(x)` converts at `precision(x)`.
2. In `IntervalArithmetic` 1.0.12, `expm1` on an `Interval{BigFloat}` of 512 bits evaluated at
   256 global bits returned a 256-bit *point* interval that excluded the true value, which an
   Arb evaluation at 1024 bits confirmed. `1 - exp(x)` on the same input stayed correct. The core
   `Seeds._defect` intersects both enclosures, so on such input the intersection was empty or
   wrong.

The certificate therefore runs `seed_defects` with the global precision set to 512 bits. It
asserts that the defect enclosures have 512 bits and positive width, and that the stored
rationals equal the exact upper ends. The core pipeline (`GapCertificates`) passes
`Interval{Float64}` to `seed_defects` and is not affected.

## Seed choice

`scripts/scan.jl` scans $(n_0, \beta_0)$ with 256-bit intervals. Float64 is useless here because
$\log Z \approx 10^4$ to $4 \times 10^5$ at the seeds while the defects are below its rounding error, and a defect
rounded to zero makes `certified_gap_bound` return a meaningless bound. With 256-bit intervals,
the best $\beta_0$ at fixed $n_0$ grows roughly linearly in $n_0$, and the bound tends to the exact
gap slowly:

| $g$ | $n_0$ (best $\beta_0$, bound) |
|---|---|
| 2 | 64 (22, 1.8511), 96 (32, 1.8942), 128 (44, 1.9188), 192 (66, 1.9432), 256 (78, 1.9498) |
| 5/4 | 96 (40, 0.4288), 128 (56, 0.4461), 144 (64, 0.4509), 192 (84, 0.4620), 256 (112, 0.4704) |
| 11/10 | 192 (88, 0.1682), 256 (120, 0.1753), 384 (176, 0.1824), 512 (240, 0.1864) |

The certified points keep the survey's $n_0$, with $\beta_0$ refined on a finer grid.
