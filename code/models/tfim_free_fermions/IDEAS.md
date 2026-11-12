# Route analysis: transverse-field Ising ring at $h = 0$ (free fermions)

Targets: the three open targets of `AGENTS.md` and `home_page/_data/records.yml` for
$`H_L = -\sum_j Z_j Z_{j+1} - g \sum_j X_j`$ on periodic rings of $L$ sites, at $g = 2$, $5/4$ and
$11/10$, with exact gap $2(g - 1)$ (paramagnetic phase). The survey estimates the seeds
$(n_0, \beta_0) = (128, 32)$, $(144, 64)$ and $(384, 161)$, far beyond exact diagonalization.

Reflection positivity is immediate: in the $\sigma^z$ basis $-h = Z \otimes Z = Z \otimes \overline{Z}$
($C = Z$) and the on-site term $o = -gX$ is real symmetric.

## Candidate routes

1. **Verified spectra (`certify_gap`, as in `tfim_longitudinal`).** Rings of $2n_0 \ge 256$ sites
   are out of reach (dense ED stops near 20 sites); small seeds ($n_0 = 6$) work but cap the
   bound well below the survey values. Rejected as the main route.
2. **Free-fermion closed form.** Jordan–Wigner maps $H_L$ to free fermions with parity-dependent
   boundary conditions (Neveu–Schwarz on even parity, Ramond on odd), so
   ```math
   Z_L(\beta) = \tfrac12\Big[\prod_{\mathrm{NS}} 2\cosh\tfrac{\beta\varepsilon_k}{2} + \prod_{\mathrm{NS}} 2\sinh\tfrac{\beta\varepsilon_k}{2}
   + \prod_{\mathrm{R}} 2\cosh\tfrac{\beta\varepsilon_k}{2} - \prod_{\mathrm{R}} 2\sinh\tfrac{\beta\varepsilon_k}{2}\Big],
   \qquad \varepsilon_k = 2\sqrt{1 + g^2 - 2g\cos k}.
   ```
   This is $O(L)$ elementary functions, so $\log Z$ at $L$ in the hundreds or thousands is
   enclosed in BigFloat interval arithmetic (MPFR, directed rounding) in milliseconds. The seed
   defects are of order $e^{-n_0/\xi}$ and $e^{-\beta_0 \Delta}$, so $\log Z \approx 10^4$–$10^5$ must be
   enclosed to far better than $10^{-30}$ absolute accuracy: Float64 cannot be used even for
   scans, 512-bit intervals are ample. Cost: linear in $n_0$, independent of $\beta_0$.
   For $g < 1$ the Ramond zero mode carries the signed energy $2(g - 1)$; all three targets have
   $g > 1$, where every $\varepsilon_k > 0$, and the code refuses $g \le 1$.
3. **Transfer-operator / expansion bounds.** Unnecessary here: the closed form is exact.

## Choice

Route 2. The closed form is classical (Lieb–Schultz–Mattis 1961, Pfeuty 1970); the tests
compare it with exact diagonalization of rings of 4 to 10 sites. The certificate reuses the core
`Seeds.seed_defects`, `certified_gap_bound` and `rational_gap_chain`.

## Outcome

Route 2 works at all three points, at the survey's seed lengths: $\Delta_L \ge 1.91897$
($g = 2$, $L \ge 128$), $0.45128$ ($g = 5/4$, $L \ge 144$) and $0.18259$ ($g = 11/10$,
$L \ge 384$). These are within $`4\text{–}10\%`$ of the exact gaps and at least the survey's
floating-point estimates. Certification takes under a second per point. Larger $n_0$ would give
more (for example $1.9498$ at $n_0 = 256$ for $g = 2$, from the scan).

Two pitfalls showed up on the way, both documented in `METHOD.md`: Julia's
`Rational{BigInt}(::BigFloat)` rounds to the global precision, and `IntervalArithmetic`'s BigFloat
`expm1` is not an enclosure when its argument carries more bits than the global precision.
