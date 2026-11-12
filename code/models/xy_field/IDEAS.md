# Route analysis: XY chain in a field with $J_x J_y < 0$ (free fermions)

Target: the open target `target-xy-field` of `home_page/_data/records.yml`,
$`H_L = \sum_j (J_x X_j X_{j+1} + J_y Y_j Y_{j+1}) - h \sum_j Z_j`$ on periodic rings of $L$ sites
with $J_x J_y < 0$. The target names no parameter point; this memo picks three.

Reflection positivity: since $\overline{Y} = -Y$,
$`-h_{\text{bond}} = (-J_x)\, X \otimes \overline{X} + J_y\, Y \otimes \overline{Y}`$, which is a positive
combination for $J_x < 0 < J_y$, with the real symmetric on-site term $o = -hZ$. The other sign
pattern $J_x > 0 > J_y$ reduces to this one by the on-site rotation $X \to Y$, $Y \to -X$, which
fixes $Z$. The certificates use $J_x = -1$.

## Candidate routes

1. **Verified spectra (`certify_gap`).** These work only for small seeds (rings of 12–20 sites),
   and the bound then falls well short of the gap. Not needed.
2. **Free-fermion closed form**, as in `code/models/tfim_free_fermions/`. Jordan–Wigner with $Z_j$ as
   the occupation parity gives hopping $\propto J_x + J_y$ and pairing $\propto J_x - J_y$, with
   ```math
   \varepsilon_k = 2\sqrt{(h - (J_x + J_y)\cos k)^2 + (J_x - J_y)^2 \sin^2 k},
   ```
   Neveu–Schwarz momenta on even parity and Ramond momenta on odd parity. The unpaired Ramond modes
   $k = 0, \pi$ carry the signed energies $2(h - (J_x + J_y))$ and $2(h + (J_x + J_y))$, and
   $`Z_L = \frac12[\prod_{\mathrm{NS}} 2\cosh + \prod_{\mathrm{NS}} 2\sinh + \prod_{\mathrm{R}} 2\cosh - \prod_{\mathrm{R}} 2\sinh]`$.
   A floating-point prototype matched dense ED at $L = 4, 6, 8$ to $10^{-14}$ for four parameter
   sets, including one in the ordered phase. In the paramagnetic phase $h > |J_x + J_y|$ every
   $\varepsilon_k$ is positive, the ground state is unique, and the gap is
   $2(h - |J_x + J_y|)$: for $J_x J_y < 0$, $|J_x + J_y| < |J_x - J_y|$, so $\varepsilon_k^2$ is concave
   in $\cos k$ and is smallest at $k = 0$ or $\pi$. The code requires this regime. Cost: linear in $L$.
3. **Ordered phase** $h < |J_x + J_y|$: the two quasi-degenerate ground states on rings rule out
   the criterion (it proves a unique ground state). Not attempted.

## Choice

Route 2 at $J_x = -1$, $h = 1$, and $J_y \in \{1/4, 1/2, 1\}$, with exact gaps $1/2$, $1$ and $2$.
The interval code mirrors `tfim_free_fermions/src/FreeFermions.jl`, including its BigFloat
precision safeguards (see that model's `METHOD.md`).

## Outcome

Route 2 works. At $n_0 = 128$ the certified bounds are $\Delta_L \ge 0.43756$ ($J_y = 1/4$),
$0.91890$ ($J_y = 1/2$) and $1.86496$ ($J_y = 1$) for every even $L \ge 128$, against exact gaps
$1/2$, $1$ and $2$. Certification takes under a second per point. The interval code mirrors
`tfim_free_fermions`, including its BigFloat precision safeguards.
