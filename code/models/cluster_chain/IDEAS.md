# Route analysis: cluster chain in a transverse field, $g = 1/2$

Target: `target-cluster-g-1-2` of `home_page/_data/records.yml`,
$`H = -\sum_j Z_{j-1} X_j Z_{j+1} - g \sum_j X_j`$ on periodic rings of $N = 2L$ sites, i.e. $L$
two-site cells, at $g = 1/2$ (cluster SPT phase, exact gap $2(1 - g) = 1$). The survey estimates
the seed $(n_0, \beta_0) = (128, 64)$ cells, float bound $0.932$.

## Reflection positivity

The bond between cells $(2c, 2c+1)$ and $(2c+2, 2c+3)$ is
$`h = -[(1 \otimes Z) \otimes (X \otimes Z) + (Z \otimes X) \otimes (Z \otimes 1)]`$, and the on-site
term is $o = -g(X \otimes 1 + 1 \otimes X)$. Neither is reflection positive in the product basis.
In the cell frame $U = (H \otimes H)\,\mathrm{CS}$, with $\mathrm{CS} = \mathrm{diag}(1, 1, 1, i)$
($H \otimes H$ has entries $\pm 1/2$, so $U$ has exact Gaussian-rational entries), the rotated
on-site term is real symmetric. The core `rp_certificate` (exact rational LDL†) then factors
the realigned bond. A screening run found that $(H \otimes H)\,\mathrm{CS}^\dagger$ works as well, while
$\mathrm{CS}(H \otimes H)$, $\mathrm{CS}^\dagger(H \otimes H)$, $\mathrm{CZ}(H \otimes H)$, $(H \otimes H)\,\mathrm{CZ}$,
$\mathrm{CS}$, $H \otimes H$ and the identity fail.

## Candidate routes

1. **Verified spectra** of cell rings: up to about 10 cells (20 sites), far below the estimated seed.
2. **Free-fermion closed form.** Jordan–Wigner with $X_j$ as the occupation parity turns
   $Z_{j-1} X_j Z_{j+1}$ into a next-nearest-neighbor Majorana coupling. The boundary terms pick
   up the fermion parity as in the transverse-field Ising chain, and the mode energies become
   $`\varepsilon_k = 2\sqrt{1 + g^2 - 2g\cos 2k}`$, with Neveu–Schwarz momenta on even parity and
   Ramond momenta on odd parity:
   $`Z_N = \frac12[\prod_{\mathrm{NS}} 2\cosh + \prod_{\mathrm{NS}} 2\sinh + \prod_{\mathrm{R}} 2\cosh - \prod_{\mathrm{R}} 2\sinh]`$.
   On even rings the Ramond modes $k = 0, \pi$ share the signed energy $2(g - 1)$, so their signs
   cancel in the $\sinh$ product, and the unsigned roots can be used for every $g \ne 1$ (the ED
   check below fixes these conventions). A floating-point prototype matched dense ED at $N = 4, 6, 8, 10$ for $g = 1/2$ and
   $g = 2$ to $10^{-14}$. The frame change $U$ does not change $Z$. Cost: linear in $N$.

## Choice

Route 2, with the interval code of `code/models/tfim_free_fermions/` adapted to $\cos 2k$ and
two-site cells.

## Outcome

Route 2 works. At the survey's seed length $n_0 = 128$ cells, $\beta_0 = 175/2$, the certified
bound is $\Delta_L \ge 0.951567$ for every even $L \ge 128$ cells (exact gap $1$, survey estimate
$0.932$). The RP factorization in the frame $(H \otimes H)\,\mathrm{CS}$ has two factors with
entries in $\{0, \pm 1, \pm i\}$ and unit weights.
