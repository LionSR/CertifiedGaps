# Method: uniform gap certificates for the Heisenberg chain in a staggered field

This document describes the reflection-positive frame, the two-site blocking that splits the seed
spectra, the certified seeds, the 20-site numerical exploration, and the history of the
certificates for the spin-1/2 Heisenberg chain in a staggered longitudinal field.

## Outcome and scope

All five target staggered fields $`h \in \{1, 3/4, 1/2, 7/20, 1/4\}`$ are certified from complete
verified 8- and 16-site ring spectra at seed length $n_0 = 8$. The alternating sublattice rotation
gives a homogeneous reflection-positive ring with a one-spin unit cell on every even ring length.
The gap bounds therefore apply to **every even chain length $L \ge 8$ sites**, establishing both
ground-state uniqueness and $\Delta_L \ge \Delta_{\mathrm{cert}}$. Each certificate is a core
certificate (`CertifiedGaps.GapCertificates`, format `CertifiedGaps ring gap certificate v1`) with
no model-local numerical code.

An exploratory numerical scan reaches $2n_0 = 20$ sites ($n_0 = 10$). Its candidate bounds are
reported below as non-certified benchmarks.

## Model, frame, and exact RP identity

We consider the periodic spin-1/2 chain
```math
H_L = \sum_{j=1}^L \mathbf{S}_j \cdot \mathbf{S}_{j+1} - h \sum_{j=1}^L (-1)^j S_j^z,
```
where the $L = 2$ ring includes both periodic bonds. Define the exact rational spin-1/2 matrices
```math
X = S^x = \begin{pmatrix} 0 & 1/2 \\ 1/2 & 0 \end{pmatrix}, \qquad
G_y = i S^y = \begin{pmatrix} 0 & 1/2 \\ -1/2 & 0 \end{pmatrix}, \qquad
Z = S^z = \begin{pmatrix} 1/2 & 0 \\ 0 & -1/2 \end{pmatrix}, \qquad U = 2G_y.
```
(In `src/common.jl`, $G_y$ is denoted `J`.) The physical bond is
$\mathbf{S} \cdot \mathbf{S} = X \otimes X - G_y \otimes G_y + Z \otimes Z$. Rotating every odd site
by the real orthogonal $U$ (a $\pi$ rotation about $y$, $UXU^\top = -X$, $UG_yU^\top = G_y$,
$UZU^\top = -Z$) turns both bond orientations into
```math
h_{\mathrm{rot}} = -(X \otimes X + G_y \otimes G_y + Z \otimes Z), \qquad
-h_{\mathrm{rot}} = \sum_{C \in \{X, G_y, Z\}} C \otimes \overline{C},
```
and the staggered field into the uniform real on-site term $o = -hZ$. The bond factors are real,
which is all that reflection positivity with $F = I$ requires; the antisymmetry of $G_y$ is
allowed. `test/runtests.jl` checks
$V_L H_L V_L^\top = \sum_j h_{\mathrm{rot}, j, j+1} + \sum_j o_j$ with
$V_L = \bigotimes_{j\ \mathrm{odd}} U_j$ exactly on rings of 2, 4 and 6 sites for every field.
The core verifies the RP identity exactly (`RingModels.verify_rp`).

## Two-site blocking of the seed spectra

The one-site rotated ring conserves only the $S^z$ parity, since its bond contains
$X \otimes X + G_y \otimes G_y = (S^+S^+ + S^-S^-)/2$. In two-site cells $(a, b)$, this term flips
both spins of a pair together: $(a, b)$ inside a cell, or $(b, a')$ across neighboring cells. The
cell charge $`q = S^z_a - S^z_b \in \{0, 1, -1, 0\}`$ is therefore an exactly conserved U(1)
charge. It equals minus the physical magnetization of the cell. The core option
`blocking = (; cell = 2, charge = [0, 1, -1, 0])` builds the exact cell model
(`RingModels.cell_model`: on-site term $h_{ab} + o_a + o_b$, bond between $b$ and the next $a$),
checks charge conservation exactly, and takes $Z_8$ and $Z_{16}$ from rings of 4 and 8 cells,
which are the same operators as the 8- and 16-site rings. Each ring splits into (cell charge,
cell momentum) blocks. Momenta $\pm k$ have conjugate blocks (the model is real), so one block
per pair is built with multiplicity 2. The certificate stores the cell model, and
`check_certificate` rebuilds it from the stored site model and compares the two.

| Ring | Blocks built | Largest block | Weighted dimension |
|---:|---:|---:|---:|
| 8 sites (4 cells) | 23 | 20 | 256 |
| 16 sites (8 cells) | 77 | 1620 | 65536 |

No on-site cell operator maps $q$ to $-q$. A global spin flip reverses the field, and the
exchange of sublattices that would restore it is a translation or reflection, not an on-site map.
The old local pipeline used spin flip composed with a one-site translation for this purpose,
so the core builds about twice as many blocks. Its cost is still 93–119 s per field.

## Enclosure and seed defects

Every block is enclosed by `CertifiedGaps.VerifiedSpectra.verified_eigen` (see
`julia/src/VerifiedSpectra.jl`). Its inputs are the exactly Hermitian `Float64` center of the
interval block matrix and an exact bound on the distance to the true block. Real blocks
(momentum $0$ and $\pi$) use the real path. The function bounds the residual and orthogonality
defect of the LAPACK eigenvectors under the IEEE binary64 + BLAS model with a priori
$\gamma_{N+8}$ or $\gamma_{2N+8}$ bounds on Frobenius norms. It then takes the smaller of the
Ostrowski and polar-factor radii and adds the matrix error by Weyl. The largest radius is
$1.2 \times 10^{-12}$ on 8 sites and $2.1 \times 10^{-8}$ on 16 sites at $h = 1$. Under the
earlier local enclosure these radii were $2.2 \times 10^{-11}$ and $2.1 \times 10^{-7}$.
`log_partition_enclosure` sums the spectra in interval arithmetic. `seed_bounds` takes the
upper ends of $1 - P_{\mathrm{s}}(8, \beta_0)$ and $1 - P_\tau(16, \beta_0)$ and rounds them up to
the grid $10^{-15}$. `certify_gap` runs the exact rational chain and the interval screen
`certified_gap_bound`.

| $h$ | $\beta_0$ | $\varepsilon_{\mathrm{s}}$ | $\varepsilon_\tau$ | Interval bound | Chain bound | 16-site ED gap | Time (s)* |
|---|---:|---|---|---:|---:|---:|---:|
| `1//1` | `21//4` | `27705127843//40000000000000` | `364207460943//500000000000000` | 1.148970072 | 1.147763329 | 1.87351193 | 93.1 |
| `3//4` | `21//4` | `585396962529//250000000000000` | `3118610756377//1000000000000000` | 0.890950432 | 0.889743689 | 1.57089206 | 113.6 |
| `1//2` | `11//2` | `1342500085157//125000000000000` | `5576810540003//500000000000000` | 0.594821280 | 0.594245333 | 1.23131976 | 118.7 |
| `7//20` | `23//4` | `510196036833//15625000000000` | `27432706400803//1000000000000000` | 0.383520437 | 0.383244983 | 0.994380622 | 113.4 |
| `1//4` | `25//4` | `2495802294179//31250000000000` | `44760635388019//1000000000000000` | 0.214491805 | 0.214238387 | 0.812031433 | 117.0 |

*Single-field times on a lightly loaded machine. The final regeneration ran under heavy
concurrent load, and the `seconds` stored in the certificates (102–519 s) reflect that load.

## Twenty-site numerical exploration and scaling

The 20-site scan was run with the retired local scan, in $(S^z, k)$ blocks under translation by
two sites. It combined full dense diagonalization for blocks of dimension $\le 2{,}500$ with
80 low-lying Ritz values from `KrylovKit` for larger blocks (up to dimension $18{,}504$, with
total weighted dimension $\sum \mathrm{mult} \times d = 1{,}048{,}576 = 2^{20}$). The current
`scripts/scan.jl` uses the core exact diagonalization in $S^z$ sectors only (Lanczos above
dimension 2500). It produces the same kind of floating-point estimate and is slower at 20 sites.
Because Krylov convergence and eigenvalue completeness are not certified, the $n_0 = 10$ rows below
are exploratory numerical estimates:

| $h$ | Best $n_0$ | $\beta_0$ | Numerical candidate bound | Numerical 20-site ED gap |
|---|---:|---:|---:|---:|
| `1//1` | 10 | 6.25 | 1.27832796 | 1.87351191 |
| `3//4` | 10 | 6.5 | 1.01075578 | 1.57089173 |
| `1//2` | 10 | 6.75 | 0.707243353 | 1.23131113 |
| `7//20` | 10 | 7 | 0.493229271 | 0.994292667 |
| `1//4` | 10 | 7.25 | 0.326952394 | 0.811528896 |
| `1//5` | 10 | 7.5 | 0.231456632 | 0.708641218 |
| `3//20` | 10 | 8 | 0.119366516 | 0.594804607 |
| `1//10` | — | — | no grid seed at $n_0 \le 10$ | 0.466432604 |

At $n_0 \le 10$, the sampled $\beta_0$ grid finds an admissible seed at $h = 3/20 = 0.15$ and none
at $h = 1/10 = 0.10$. The table below records $\max_\beta \min(P_{\mathrm{s}}(n_0, \beta), P_\tau(2n_0, \beta))$
and the ground-state spatial energy defect $\delta(n_0) = E_0(2n_0) - 2E_0(n_0)$ across all scanned
fields and seed lengths $`n_0 \in \{2, 4, 6, 8, 10\}`$:

| $h$ | $n_0$ | $\max_\beta \min(P_{\mathrm{s}}, P_\tau)$ | $\delta(n_0) = E_0(2n_0) - 2E_0(n_0)$ |
|---|---:|---:|---:|
| `1//1` | 2 | 0.529873172 | 0.505022849 |
| `1//1` | 4 | 0.94067001 | 0.0221827722 |
| `1//1` | 6 | 0.99334286 | 0.00159103285 |
| `1//1` | 8 | 0.999275644 | 0.00013128857 |
| `1//1` | 10 | 0.999923253 | 1.14962961e-05 |
| `3//4` | 2 | 0.447640968 | 0.613999064 |
| `3//4` | 4 | 0.890324463 | 0.0389149028 |
| `3//4` | 6 | 0.983514815 | 0.00385786603 |
| `3//4` | 8 | 0.997555655 | 0.000440747681 |
| `3//4` | 10 | 0.999636029 | 5.36822693e-05 |
| `1//2` | 2 | 0.399066156 | 0.754873673 |
| `1//2` | 4 | 0.800725359 | 0.0769369751 |
| `1//2` | 6 | 0.952122862 | 0.0113905614 |
| `1//2` | 8 | 0.988846519 | 0.00192196739 |
| `1//2` | 10 | 0.997543855 | 0.000347058281 |
| `7//20` | 2 | 0.361338696 | 0.850802846 |
| `7//20` | 4 | 0.699446862 | 0.124633121 |
| `7//20` | 6 | 0.893328792 | 0.0251937055 |
| `7//20` | 8 | 0.967347486 | 0.00566434475 |
| `7//20` | 10 | 0.990324488 | 0.00135842657 |
| `1//4` | 2 | 0.31998312 | 0.913017541 |
| `1//4` | 4 | 0.602523494 | 0.177017263 |
| `1//4` | 6 | 0.818031626 | 0.0466268272 |
| `1//4` | 8 | 0.926399235 | 0.0132591633 |
| `1//4` | 10 | 0.971872848 | 0.00397943173 |
| `1//5` | 2 | 0.303223146 | 0.940972384 |
| `1//5` | 4 | 0.560054301 | 0.2117524 |
| `1//5` | 6 | 0.757748899 | 0.0653997426 |
| `1//5` | 8 | 0.881241013 | 0.0214187006 |
| `1//5` | 10 | 0.947706946 | 0.00732718927 |
| `3//20` | 2 | 0.289667018 | 0.965093873 |
| `3//20` | 4 | 0.495941276 | 0.252166166 |
| `3//20` | 6 | 0.674748062 | 0.0934440211 |
| `3//20` | 8 | 0.814781105 | 0.0361784563 |
| `3//20` | 10 | 0.902221022 | 0.0144183736 |
| `1//10` | 2 | 0.279689865 | 0.983874689 |
| `1//10` | 4 | 0.442024104 | 0.295440509 |
| `1//10` | 6 | 0.569803257 | 0.134401906 |
| `1//10` | 8 | 0.7028359 | 0.0639327614 |
| `1//10` | 10 | 0.80680422 | 0.0308278889 |

### Asymptotic $h^{2/3}$ scaling comparison

In the continuum sine-Gordon description at weak staggered field $h \ll 1$, the spectral gap and
inverse correlation length scale as $\Delta \propto h^{2/3}$ (up to logarithmic corrections).
Calibrating to the $h = 1$ seed $(n_0, \beta_0) = (8, 5.25)$ gives the constant-defect-quality
scaling locus $n_0(h) \approx 8 h^{-2/3}$ and $\beta_0(h) \approx 5.25 h^{-2/3}$:

| $h$ | Scaled $n_0 \approx 8 h^{-2/3}$ | Next even $n_0$ | Scaled $\beta_0 \approx 5.25 h^{-2/3}$ | Certified $n_0$ |
|---|---:|---:|---:|---:|
| `1//1` | 8 | 8 | 5.25 | 8 |
| `3//4` | 9.69130983 | 10 | 6.35992207 | 8 |
| `1//2` | 12.6992084 | 14 | 8.33385552 | 8 |
| `7//20` | 16.1081114 | 18 | 10.5709481 | 8 |
| `1//4` | 20.1587368 | 22 | 13.229171 | 8 |

Because the admissible region allows seed defects up to $\approx 0.12731$ (whereas the $h = 1$ seed
has defects $\sim 7 \times 10^{-4}$), $n_0 = 8$ remains sufficient all the way down to $h = 1/4$,
with the certified bound-to-ED-gap ratio decreasing from $0.613$ at $h = 1$ to $0.264$ at $h = 1/4$.
At $h = 1/10$, fitting $\log \delta(n_0) \approx 0.19944 - 0.36810 n_0$ over $`n_0 \in \{6, 8, 10\}`$
gives a decay scale of $2.72$ sites and predicts a first admissible seed near
$(n_0, \beta_0) \approx (14, 17.0)$ (requiring a 28-site doubled ring).


## History of the certificates

The five certificates were first produced by a model-local pipeline in the physical frame. It
built $(S^z, k)$ blocks under translation by two sites, using spin flip composed with a one-site
translation to identify $p \leftrightarrow L - p$. Its own residual and polar-factor enclosure used
conservative $\gamma_{32d+32}$ constants and $`\|A\|_2 \le L(3/4 + h/2)`$. It had its own $\log Z$
and a fresh spectral replay and direct interval audit of small blocks. An independent
shared-core certificate at $h = 1$ in the one-site frame used only the $S^z$ parity
(about 12 min, $\Delta_L > 1.147699$). The local pipeline was retired in favor of the core.
The regenerated certificates have smaller defects at every point because the core radii are
about ten times smaller:

| $h$ | Old $\varepsilon_{\mathrm{s}}$ | New $\varepsilon_{\mathrm{s}}$ | Old $\varepsilon_\tau$ | New $\varepsilon_\tau$ | Old interval bound | New interval bound |
|---|---:|---:|---:|---:|---:|---:|
| $1$ | $6.93687 \times 10^{-4}$ | $6.92628 \times 10^{-4}$ | $7.32648 \times 10^{-4}$ | $7.28415 \times 10^{-4}$ | 1.148261490 | 1.148970072 |
| $3/4$ | $2.34252 \times 10^{-3}$ | $2.34159 \times 10^{-3}$ | $3.12234 \times 10^{-3}$ | $3.11861 \times 10^{-3}$ | 0.890782204 | 0.890950432 |
| $1/2$ | $1.074085 \times 10^{-2}$ | $1.074000 \times 10^{-2}$ | $1.115701 \times 10^{-2}$ | $1.115362 \times 10^{-2}$ | 0.594785233 | 0.594821280 |
| $7/20$ | $3.265334 \times 10^{-2}$ | $3.265255 \times 10^{-2}$ | $2.743591 \times 10^{-2}$ | $2.743271 \times 10^{-2}$ | 0.383508307 | 0.383520437 |
| $1/4$ | $7.986645 \times 10^{-2}$ | $7.986567 \times 10^{-2}$ | $4.476387 \times 10^{-2}$ | $4.476064 \times 10^{-2}$ | 0.214486106 | 0.214491805 |
