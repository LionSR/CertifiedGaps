# Spin-2 Heisenberg chain: thermal-MPO investigation

This Julia/ITensors project investigates the periodic Hamiltonian

\[
H_L=\sum_{j=1}^L \mathbf S_j\cdot\mathbf S_{j+1},\qquad
\mathbf S^2=6I,\quad J=1,\quad \mathbf S_{L+1}=\mathbf S_1.
\]

It implements the spin-2 transfer algebra, periodic DMRG, and a numerical
two-purity scan. **A uniform spin-2 gap has not been certified.** The analytic
transfer construction extends to spin two; the missing ingredients are rigorous
finite-scale purity estimates. The runs here establish numerical baselines and
check the implementation against independent small-system calculations.

## Run

The project uses Julia 1.12, ITensors 0.9.32, and ITensorMPS 0.4.2. Exact resolved
dependencies are in `Manifest.toml`. All commands below run from the repository root.

```sh
julia --startup-file=no --project=models/spin2 -e 'using Pkg; Pkg.instantiate()'
julia --startup-file=no --project=models/spin2 models/spin2/test/runtests.jl
julia --startup-file=no --project=models/spin2 models/spin2/scripts/run.jl \
  models/spin2/configs/dmrg-pilot.toml models/spin2/results/dmrg-pilot.json
julia --startup-file=no --project=models/spin2 models/spin2/scripts/run.jl \
  models/spin2/configs/thermal-pilot.toml models/spin2/results/thermal-pilot.json
```

On the initial machine, Julia is installed outside `PATH` at
`/Users/siruilu/.julia/juliaup/julia-1.12.4+0.aarch64.apple.darwin14/Julia-1.12.app/Contents/Resources/julia/bin/julia`.
Use that executable in place of `julia` if needed.

Format only the Julia files you changed (the settings are in `.JuliaFormatter.toml`:
margin 100, `for … in`). Never format the whole tree: several results record SHA-256
hashes of the sources that produced them.

```sh
julia --startup-file=no --project=models/spin2 -e 'using JuliaFormatter; format_file("models/spin2/src/exact.jl")'
```

## Code layout

`spin2` is the package `Spin2Gap`; it depends on the shared core package
`CertifiedGaps` in `../julia` (a path dependency recorded in `Manifest.toml`, so
`Pkg.instantiate()` is enough). Scripts and tests load code with `using Spin2Gap`
and `import CertifiedGaps`; no script includes a `src/` file. Scripts still include
other scripts where one is a library of another, and are still run from the repository
root with `--project=models/spin2`.

| module | file(s) | content |
|---|---|---|
| `Spin2Gap` | `src/Spin2Gap.jl`, `model.jl`, `exact.jl`, `dmrg.jl`, `thermal.jl` | S=2 site type, `bonds`, `hamiltonian`, exact spectra, DMRG, thermal curves, `provenance`, `write_result` |
| `Spin2Gap.Trials` | `src/Trials.jl`, `src/trial_certificate.jl` | SU(2) trial templates and exact integer energy certificates (was `Spin2Trials`) |
| `Spin2Gap.Variational` | `src/Variational.jl` | reduced variational energies and optimizers (was `Spin2Variational`) |
| `Spin2Gap.Certificate` | `src/Certificate.jl` | reduced exact/ball trial certificates (was `Spin2Certificate`) |
| `Spin2Gap.ThermalBounds` | `src/ThermalBounds.jl` | rational Trotter-layer thermal bounds (was `Spin2ThermalBounds`) |
| `Spin2Gap.PolynomialMPO` | `src/PolynomialMPO.jl` | spin-2 polynomial and quartic SOS seeds (was `Spin2PolynomialMPO`) |

The certified MPO toolkit moved verbatim into the core, `CertifiedGaps.CertifiedMPO`;
each submodule keeps the old aliases (`F`, `MN`, `MS`, `CB`, `PG`):

| former spin2 module | now |
|---|---|
| `Spin2MPONorms` | `CertifiedGaps.CertifiedMPO.Norms` |
| `Spin2MPOSquaring` | `CertifiedGaps.CertifiedMPO.Squaring` |
| `Spin2MPOProducts` | `CertifiedGaps.CertifiedMPO.Products` |
| `Spin2MPOOperatorError` | `CertifiedGaps.CertifiedMPO.OperatorError` |
| `Spin2MPOGaugeSearch` | `CertifiedGaps.CertifiedMPO.GaugeSearch` |
| `Spin2AdjointTrace` | `CertifiedGaps.CertifiedMPO.AdjointTrace` |
| `Spin2SchattenBridge` | `CertifiedGaps.CertifiedMPO.SchattenBridge` |
| `Spin2PositiveBounds` | `CertifiedGaps.CertifiedMPO.PositiveBounds` |
| `Spin2PolynomialBounds` | `CertifiedGaps.CertifiedMPO.PolynomialBounds` |
| generic part of `Spin2PolynomialMPO` | `CertifiedGaps.CertifiedMPO.MPOAlgebra` (re-bound in `Spin2Gap.PolynomialMPO`) |

Other shared code now comes from the core: `spin_matrices()` and `real_generators()` are
`Spins.spin_matrices(2)` and `Spins.cartesian_generators(2)`; `Trials.ladders` and
`Trials.WEIGHTS` are `Spins.integer_ladders` and `Spins.binomial_weights(2)`;
`sector_hamiltonian`, `exact_spectrum` and `log_partition` use
`ExactDiagonalization` (charge 2m per site, spin-flip doubling); `thermal_curve` wraps
`Purification.purified_log_partition` (its β = 0 row now also has `"seconds" => 0.0`);
`provenance` and `write_result` use `Results`; `scripts/screen.jl` uses
`Recurrence.beta_window`. Pretty-printed JSON results are written atomically by
`write_result`, which adds a trailing newline; compact witness files are written as before.

`scripts/bootstrap_screen.jl` and `scripts/one_trace_target.jl` no longer include
`models/spin1/experiments/optimize_bounds.jl`: the exact scalar chain is
`CertifiedGaps.Recurrence.rational_gap_chain` (a verbatim port of spin1's `gap_chain`,
grid 10^512) with `ExactArithmetic.frac`/`first_true`. Their outputs are unchanged except
that `scalar_backend_sha256` in `one-trace-targets.json` now hashes
`julia/src/recurrence/rational.jl` (stored files from before the move hash the spin1 file).
Hash fields that scripts record for modules (`norm_module_sha256`,
`squaring_module_sha256`, `bounds_source_sha256`, `source_sha256` lists) likewise point at
the new files.

The tests of the moved modules now live in the core:
`julia/test/certified_mpo/{mpo_norms,mpo_squaring,mpo_products,charge_blocks,mpo_operator_error,mpo_gauge_search,adjoint_trace,schatten_bridge,positive_bounds,polynomial_bounds}.jl`
(each runs standalone with `--project=julia`, or all with
`julia --project=julia -e 'import Pkg; Pkg.test()'`). The spin2 tests are `runtests.jl`
plus the standalone files `trials.jl`, `variational.jl`, `certificate.jl`,
`thermal_bounds.jl`, `polynomial_mpo.jl`, `l2_bounds.jl`, `mpo_adjoint_flow.jl` and
`one_trace.jl`, each run with `--project=models/spin2`.

To reproduce the tighter runs, use `configs/dmrg-refinement.toml` and
`configs/thermal-refinement.toml` with `scripts/run.jl`, writing respectively to
`results/dmrg-refinement.json` and `results/thermal-refinement.json` (all paths
prefixed with `spin2/` when invoked from the repository root).

The runner writes JSON after each completed parameter combination, with a
`complete` field distinguishing finished and interrupted runs. Each result records
the configuration, dependency versions, thread counts, source hashes, sweep
histories or thermal normalization data, and `certified=false`. Rerunning with the
same output path replaces it. Give separate runs separate output filenames.

## Why the exact transfer argument extends

Real symmetric traceless 3-by-3 matrices form a five-dimensional spin-2 space.
For the real antisymmetric spin-1 rotation generators \(J_\alpha\), let

\[
G_\alpha Q=J_\alpha Q-QJ_\alpha,\qquad S^\alpha=-iG_\alpha.
\]

In an orthonormal real basis of this space, the \(G_\alpha\) are real
antisymmetric matrices. They satisfy
\([G_x,G_y]=-G_z\), \(-\sum_\alpha G_\alpha^2=6I\), and
\(\|G_\alpha\|=2\). `real_generators()` builds this representation independently
of the conventional magnetic-basis operators used by ITensors.

Set \(A_\alpha=G_\alpha/2\). Then \(\|A_\alpha\|=1\) and

\[
-\mathbf S_j\cdot\mathbf S_{j+1}
=4\sum_\alpha A_{\alpha,j}A_{\alpha,j+1}.
\]

Use the same ordered bond histories as in the spin-1 note, now with weight
\(4^d\) for a history of \(d\) events. The history measure has total mass
\(e^{12\beta}\). The local trace kernel is real, symmetric under exchange of
its incident histories, and bounded by 5. Consequently its integral operator
\(X_\beta\) is self-adjoint and Hilbert--Schmidt. The absolutely convergent bond
expansion yields

\[
Z_L(\beta)=\operatorname{Tr}X_\beta^L.
\]

For even \(L\), \(X_\beta^L/Z_L(\beta)\) is positive with unit trace. The
existing spatial/thermal squaring, coupling, and interpolation argument therefore
applies without a change in its scalar constants. In particular, **if** an even
\(n\ge4\), \(\beta>0\), and \(0<\varepsilon\le0.01\) satisfy the certified bounds

\[
P_s(n,\beta)=\frac{Z_{2n}(\beta)}{Z_n(\beta)^2}\ge1-\varepsilon,
\qquad
P_\tau(2n,\beta)=\frac{Z_{2n}(2\beta)}{Z_{2n}(\beta)^2}\ge1-\varepsilon,
\]

then every even \(L\ge n\) has a unique ground state and

\[
\Delta_L>\frac{\log[1/(5\varepsilon)]}{\beta}.
\]

This is a conditional specialization of the proof in
[`thermal-mpo-gap-note.tex`](../outputs/thermal-mpo-gap-note.tex), not a claim that
spin-2 seed inequalities have been established.

## What the numerical calculations measure

### DMRG

`dmrg_levels` uses a custom five-state `S=2` site with conserved total \(S^z\).
ITensor's QN value is \(2S^z\). The Hamiltonian includes the closing \((L,1)\)
bond. An open-boundary MPS representation is used for this **periodic physical
Hamiltonian**; no claim about a periodic-MPS algorithm is intended.

For each chain it calculates the lowest state in \(S^z=0\), the lowest state in
\(S^z=1\), and optionally an excited \(S^z=0\) state orthogonal to the first.
The latter avoids assuming the lowest excitation must be a triplet. It reports
the physical expectation of \(H\), separately from the penalty objective,
energy variance, final energy change, discarded weights, and ground-state overlap.
Small residuals establish proximity to some eigenvalue, not completeness of the
low-energy spectrum. Differences of approximate variational energies are not
rigorous lower bounds on a gap.

The pilot compares independent runs at bond dimensions 64 and 128 and lengths
4, 8, and 16. It saves MPS checkpoints under `work/spin2/checkpoints/` for later
analysis or variational-certificate work. These Julia serialization files require
the pinned environment and should only be loaded from trusted local runs.

### Thermal MPO

`thermal_curve` evolves an MPO amplitude
\(A(\beta)=e^{-\beta H/2}\), starting from \(I/\sqrt{5^L}\). A forward/reverse
bond sweep applies half-step gates \(e^{-\delta\beta h/4}\), giving a symmetric
second-order product formula. It evolves the ket legs only and includes the
periodic closing bond. The norm is the Hilbert--Schmidt norm, so accumulated
normalizations recover

\[
\log Z_L(\beta)=L\log 5+2\sum_{\text{steps}}\log(\text{norm growth}).
\]

The trace includes **all magnetization sectors**. QN conservation is used only
to make tensors sparse; there is no restriction of the thermal trace to \(S^z=0\).
Ratios are evaluated in logarithmic form, with the thermal seed taken at length
\(2n\). Values outside the physical purity range are flagged, not clipped.

The thermal pilot crosses two time steps with two bond dimensions, rather than
changing both at once. It benchmarks length 4 against full diagonalization and
uses lengths 4 and 8 to exercise the entire seed calculation. These short,
high-temperature runs are implementation checks, not likely certification scales.
Trotter error and MPO truncation both require rigorous bounds before numerical
partition functions can become certificates.

### Independent checks

`sector_hamiltonian` constructs matrices directly in the magnetic product basis,
without using ITensor's MPO. `exact_spectrum` diagonalizes every magnetization
sector (using spin-flip symmetry) and includes all multiplicities. It is limited
to six sites. Here “exact diagonalization” means the full untruncated Hilbert
space; eigenvalues are still ordinary Float64 approximations.

Tests check spin normalization, the real representation, the two-spin spectrum,
open and periodic MPO spectra, the analytically known four-site values
\(E_0=-20\), \(E_1-E_0=1\), DMRG residuals and orthogonality, full thermal
normalization, and second-order time-step convergence.

## Next mathematical and computational milestones

An inexpensive necessary-condition screen helps reject unsuitable starting
lengths before thermal evolution. Define
\(\delta E_n=E_0(2n)-2E_0(n)\). If both purities are at least \(p>1/2\), then

\[
\beta\Delta_{2n}\ge\log\frac{p}{1-p},\qquad
\beta\delta E_n\le-2\log p.
\]

Indeed, thermal purity at least \(p\) gives ground weight \(w_0\ge p\), hence
\(e^{-\beta\Delta_{2n}}\le(1-p)/p\). Writing
\(Z_L=e^{-\beta E_0(L)}R_L\), one has \(R_n\ge1\) and \(R_{2n}\le1/p\).
Thus \(p\le P_s=e^{-\beta\delta E_n}R_{2n}/R_n^2\le e^{-\beta\delta E_n}/p\).
For positive \(\delta E_n\), the necessary window is

\[
\frac{\log[p/(1-p)]}{\Delta_{2n}}
\le\beta\le\frac{-2\log p}{\delta E_n}.
\]

At \(p=0.99\), it requires \(\delta E_n/\Delta_{2n}\lesssim0.004375\).
`scripts/screen.jl` evaluates this window from the highest-bond-dimension DMRG
rows available at paired lengths. With approximate energies it remains a
heuristic; rigorous rejection would require suitable energy/gap enclosures.
Passing this necessary screen does not establish either purity inequality.

```sh
julia --startup-file=no --project=models/spin2 models/spin2/scripts/screen.jl \
  models/spin2/results/seed-screen.json \
  models/spin2/results/dmrg-pilot.json models/spin2/results/dmrg-refinement.json
```

The next steps are:

1. Converge periodic low-energy calculations at increasing lengths and bond
   dimensions, with repeated starts where needed. Do not extrapolate a bulk gap
   from short chains.
2. Search \((n,\beta)\) using **both** purities. Cooling improves thermal purity
   but does not automatically improve spatial purity at fixed length.
3. Replace numerical seed estimates by certified inequalities. The original
   spin-1 proof uses finite thermal moments and variational trial states; a
   spin-2 adaptation must rebuild these certificates and their error bounds.
4. Apply the scalar doubling theorem only once both rigorous seeds pass.

For seed intervals \(\ell_L(\beta)\le\log Z_L(\beta)\le u_L(\beta)\), the
required lower bounds follow from

\[
\log P_s\ge\ell_{2n}(\beta)-2u_n(\beta),\qquad
\log P_\tau\ge\ell_{2n}(2\beta)-2u_{2n}(\beta).
\]

This is the concrete interface a later certified thermal computation must supply.
DMRG gap estimates alone cannot supply it.

## Initial results (October 7, 2026)

The validation suite passes all 36 checks. Raw numerical data, configurations,
source hashes, and convergence diagnostics are saved in
[`dmrg-pilot.json`](results/dmrg-pilot.json),
[`dmrg-refinement.json`](results/dmrg-refinement.json),
[`thermal-pilot.json`](results/thermal-pilot.json), and
[`thermal-refinement.json`](results/thermal-refinement.json).

The periodic low-energy estimates at the largest tested bond dimensions are:

| Length | Bond dimension | Ground energy | Gap estimate | Ground residual norm |
| --- | --- | --- | --- | --- |
| 4 | 128 | -20.00000000 | 1.00000000 | At floating-point roundoff |
| 8 | 512 | -38.52186914 | 0.54556928 | 0.000044 |
| 16 | 512 | -76.37583432 | Approximately 0.3109 | 0.0265 |

The eight-site excitation channels agree to about \(3\times10^{-12}\), while
raising its bond dimension from 256 to 512 changes the gap estimates by about
\(1.8\times10^{-7}\). These are convergence observations, not certified error
bars. The sixteen-site estimate is preliminary: small sweep-to-sweep energy
changes coexist with appreciable residuals. None of these short-chain gaps is
an estimate of the infinite-chain gap.

For the thermal calculation at base length \(n=4\), bond dimension 128, and
\(\delta\beta=0.025\):

| Inverse temperature | Spatial purity at n=4 | Thermal purity at 2n=8 |
| --- | --- | --- |
| 0.1 | 0.9911683 | 0.0000063016 |
| 0.2 | 0.9001105 | 0.0000467644 |

No tested point satisfies both 0.99 thresholds. The largest absolute error in
the four-site log partition function over the sampled temperatures is about
\(2.8\times10^{-4}\) at the finer time step. Increasing bond dimension from 32
to 128 reduces this error substantially, but time-step refinement alone is not
monotone because repeated truncation also contributes. These runs should not be
treated as interval enclosures.

The necessary-condition screen explains the difficulty of the smallest base
length: at \(n=4\), the refined energies require approximately
\(\beta\ge8.42\) for thermal purity but \(\beta\le0.0136\) for spatial purity.
The intervals are incompatible. See [`seed-screen.json`](results/seed-screen.json)
for both paired-length screens. This supports moving to larger lengths before
attempting low-temperature thermal computations; it does not rule out the
spin-2 gap or certify any spectral bound.

## Exact variational inputs and improved bootstrap targets (October 8, 2026)

The next implementation now searches SU(2)-structured periodic MPS directly in Julia. `src/Trials.jl` constructs rank-two tensor intertwiners from exact highest-weight equations. The virtual space can contain repeated integer-spin multiplets. A small transfer matrix on their multiplicity spaces makes the numerical energy search inexpensive. Numerical optimization is only a candidate generator; it does not certify a ground energy or a gap by itself.

The current certified candidate uses `V₀ ⊕ 2V₁ ⊕ 2V₂ ⊕ V₃`, bond dimension **24**, with **29 integer reduced coefficients** after quantization. Its physical tensors are `A_m = sqrt(binomial(4,2+m)) P_m` with integer matrices `P_m`. These square roots cancel from the double-layer norm and energy insertions. `src/trial_certificate.jl` checks an explicit overflow bound before native integer construction, then uses Nemo/FLINT arbitrary-precision integer powers.

Fresh exact contractions establish **trial energy per bond < −4.7595 at N = 128, 256, 512**. Therefore the shifted partition function with `a=4.7595` exceeds one at those lengths for every positive inverse temperature. This is a finite variational certificate, **not a spin-2 gap certificate**. The 512-site trial energy is approximately −4.7595993033. The numerical bulk reference −4.761249(6) indicates substantial energy-accuracy work remains.

```sh
# Candidate search, optional when using the saved integer witness:
julia --startup-file=no --project=models/spin2 models/spin2/scripts/search_trial.jl
julia --startup-file=no --project=models/spin2 models/spin2/scripts/refine_trial.jl
julia --startup-file=no --project=models/spin2 models/spin2/scripts/polish_trial.jl

# Exact replay of the saved integer witness, without numerical optimization:
julia --startup-file=no --project=models/spin2 models/spin2/scripts/certify_trial.jl \
  models/spin2/results/trial-replay.json models/spin2/results/trial-certificate.json
julia --startup-file=no --project=models/spin2 models/spin2/test/trials.jl

# Model-independent scalar targets plus heuristic screens of the saved DMRG data:
julia --startup-file=no --project=models/spin2 models/spin2/scripts/bootstrap_screen.jl
```

The certificate stores integer parameters so the tests reconstruct the state without trusting floating-point search coefficients. A direct four-site integer calculation agrees with the transfer contraction, and an ITensors physical MPO/MPS energy agrees within 7.2e-15. All **47 new checks pass**. The exact run for all three long lengths, including construction and the ITensors control, took about 67 seconds on the current machine; this excludes candidate optimization and Julia startup.

The improved coupled bootstrap also admits targets such as spatial purity 0.95 and temporal purity 0.83. The corresponding necessary `deltaE/gap` screen is about 0.14986, versus 0.00437435 under the older 0.99/0.99 target. The existing length-4 and length-8 base scales fail even this relaxed screen. **The missing step is still certified thermal purity inputs.** Numerical gap differences and variational energy witnesses cannot replace those inputs.

The search output shows why virtual multiplicities matter: one copy each of virtual spins 0–4 gives energy near −4.75026, while repeated low-spin multiplets reach about −4.75960 at D=24. These searches were not converged to a demonstrated optimum; saved gradients and results should be treated as exploratory. Outputs are in `results/trial-search.json`, `trial-refinement.json`, `trial-polished.json`, `trial-certificate.json`, and `bootstrap-screen.json`.

## Analytic optimization and reduced certificates (October 8, continued)

`src/Variational.jl` differentiates the energy including both transfer fixed
points through a reduced resolvent. The energy channel is factored into small
multiplicity matrices, avoiding the expanded quartic polynomial. Independent
full physical contractions and five-point finite differences check the reduced
energy and derivatives. The optimizer remains a numerical candidate generator.

Completed searches give:

| Virtual bond dimension | Numerical energy per bond |
| --- | --- |
| 40 | -4.7608215778 |
| 64 | -4.7610606436 |
| 89 | -4.7611948539 |
| 137 | -4.7612284940 |
| 197 | -4.7612356898 |

A subsequent search using a preconditioner built from the transfer fixed points
improves D=197 to **-4.7612412068**, with an analytic gradient of about 2.3e-7.
Its coefficients and batch history are in `trial-preconditioned.json`; this
candidate has not gone through the exact finite-energy verifier.

The subsequent D=295 search reaches **-4.7612463445** in
`results/trial-295-refined.json`. The search now restricts the real transfer
fixed-point calculation to symmetric multiplicity matrices and uses linear
solves for its derivative. A common preconditioner floor handles nearly unused
virtual multiplets. Independent full physical energies and derivatives agree;
all 108 variational tests pass. Numerical optimization remains separate from
the finite-energy certificate.

These are the saved values in `trial-analytic.json` and `trial-large.json`, not
claims of optimization convergence. A separate earlier survey in
`../survey/models/spin2_heisenberg/` reports still better iDMRG energies; those
results do not provide the integer witness or rigorous thermal bounds needed here.

The new verifier `src/Certificate.jl` reconstructs the trial from integer
coefficients, checks exact tensor covariance and complete irreducible bases,
checks every highest-weight image, and contracts the resulting rational matrices
using outward-rounded Arb arithmetic. It uses native integer arithmetic only
after proving an overflow bound, with an arbitrary-precision fallback. At D=137,
the largest reduced matrix is 529 by 529 instead of a full 18769 by 18769 matrix.
The first run took about 489 seconds including exact reduction and all three
long-chain contractions.

`results/trial-137-certificate.json` certifies the following strict intervals
for the **finite periodic trial energy per bond**:

| Length | Lower endpoint | Upper endpoint |
| --- | --- | --- |
| 128 | -4.761268026389 | -4.761268026385 |
| 256 | -4.761228827807 | -4.761228827803 |
| 512 | -4.761228489623 | -4.761228489619 |

In particular, all three trial energies are below **-4.76122 per bond**.
This improves the previous certified shift from 4.7595 to 4.76122, while
remaining a finite variational statement. It does not certify a gap or an
infinite-volume energy enclosure.

```sh
julia --startup-file=no --project=models/spin2 models/spin2/test/variational.jl
julia --startup-file=no --project=models/spin2 models/spin2/test/certificate.jl
julia --startup-file=no --project=models/spin2 models/spin2/scripts/large_trial.jl
julia --startup-file=no --project=models/spin2 models/spin2/scripts/certify_reduced_trial.jl \
  models/spin2/results/trial-137-replay.json models/spin2/results/trial-137-certificate.json 137 replay
```

### A rigorous thermal upper bound, currently at reference sizes

`src/ThermalBounds.jl` implements a different thermal route. Write the even
periodic Hamiltonian as `H=A+B`, with disjoint odd and even bond layers. The
Araki--Lieb--Thirring inequality and the product-formula limit give

\[
Z_N(\beta)\le
\operatorname{Tr}\left(e^{-\beta A/m}e^{-\beta B/m}\right)^m
\quad(m\in\mathbb N).
\]

The right side decreases when the step count is multiplied by an integer.
Thus a rigorous contraction of this product formula gives an **upper bound**
without a separate Trotter-error estimate. This statement follows from the
[Araki--Lieb--Thirring trace inequality](https://arxiv.org/abs/math/0701129);
it does not make an ordinary truncated MPO contraction rigorous.

The current reference implementation uses exact rational two-site spectral
projectors, all magnetization sectors, and Arb matrix products. Its 42 checks
include comparison with the independently known complete four-site spectrum
and strict monotonicity of the computed upper bounds. It deliberately limits
full-sector calculations to four and six sites. The returned ball encloses
the **Trotter upper bound**, not the true partition function on both sides.

```sh
julia --startup-file=no --project=models/spin2 models/spin2/test/thermal_bounds.jl
```

The outstanding proof step is a rigorous contraction of the thermal upper
bound at useful spin-2 lengths and temperatures. Tensor compression must carry
a proven error bound; a small discarded weight or agreement between two bond
dimensions alone is insufficient. No uniform spin-2 gap has been certified.

### Entrywise positive compression experiment

The Marshall transformation makes the thermal gates nonnegative entrywise.
This permits a different error certificate: replace a nonnegative matrix `M`
by nonnegative rational factors `U,V` satisfying `M <= U*V` entrywise. Every
subsequent contraction with nonnegative tensors preserves the upper bound.
`CertifiedGaps.CertifiedMPO.PositiveBounds` (`julia/src/CertifiedMPO/PositiveBounds.jl`) constructs candidate factors numerically, quantizes
them, adds a rank-one residual envelope, and verifies the inequality exactly.
Its 91 checks also exercise contraction monotonicity and nonuniform physical
metrics. Candidate convergence is not an acceptance condition.

`scripts/positive_compression_probe.jl` checks this idea on a four-site thermal
amplitude, resolved into charge blocks. At total bond dimension 217, the first
unweighted majorants inflated the Trotter trace by factors 6.30 at beta=1 and
4.50 at beta=4. Optimizing with the physical metric and penalizing underestimation
reduces these factors to **1.0966** and **1.3079**, respectively. The exact
majorant checks pass, but the bounds remain too loose for a useful repeated
compression scheme. These are four-site reference experiments, not large-chain
purity certificates. Raw results are in `positive-compression-probe.json` and
`positive-compression-weighted.json`.

Increasing the cap to 64 per charge block gives total bond dimension 569
(the full cut has dimension 625) and factors 1.00899 at beta=1 and 1.00844
at beta=4, saved in `positive-compression-highrank.json`. This is tighter
but still offers only a small compression at this reference size.

The repeated test in `results/repeated-positive-compression.json` exposes a
larger loss: five squarings from beta=1/8 to beta=4, compressing at all six
stages, inflate the trace by **4.9194** at D=367 and **3.2550** at D=569.
Adaptive residual envelopes give only modest improvements in a single
compression. A covering LP initially suffered from badly scaled constraints;
equilibration and retaining the better repaired candidate prevent the observed
large regressions. These variants have not supplied useful long-chain bounds.

### Certified norm and spectral compression probes

Two new four-site probes allow signed SVD factors while certifying their error.
`scripts/repeated_l2_probe.jl` preserves charge and Hermiticity exactly, using
transpose-parity blocks so that Hermitian projection does not double the rank.
It rounds coefficients to rationals and checks the full physical
Hilbert--Schmidt residual. If the current error is at most `epsilon`, the
current approximate operator norm is at most `b`, and the next compression
and rounding error is at most `delta`, the next error is at most
`delta + epsilon*(2b+epsilon)`. This follows by expanding `(B+E)^2`.

`scripts/loewner_compression_probe.jl` instead maintains domination of the
ordered eigenvalues of the exact positive amplitude. If
`norm(M^2-C) <= delta`, the replacement `C+delta*I` dominates `M^2` in the
Loewner order. Min--max and positivity propagate eigenvalue domination, and
the final squared Hilbert--Schmidt norm is a thermal upper bound. This argument
does not assume that operator squaring preserves the Loewner order.

At beta=4 after the same five squarings, the completed results are:

| Method | Represented rank bound | Certified upper / uncompressed Trotter trace |
| --- | --- | --- |
| Entrywise majorant | 367 | 4.91942725 |
| Hilbert--Schmidt error | 358 | 1.23237600 |
| Spectral majorant | 359 | 1.20390020 |
| Entrywise majorant | 569 | 3.25495415 |
| Hilbert--Schmidt error | 564 | 1.00841069 |
| Spectral majorant | 565 | 1.00388214 |

The full cut has dimension 625. The new data are in
`results/repeated-l2-symmetric.json` and `results/loewner-compression.json`.
The 14 norm tests compare the physical metric, Hermitian decomposition,
operator-norm witness, interval squaring, and perturbation bound with
independent small matrices.

The operator-norm witness currently enumerates the four-site sectors and
verifies a weighted absolute-row-sum bound using strictly positive rational
vectors. **A sufficiently sharp scalable counterpart for a long-chain MPO is still missing.**
These experiments improve the reference compression bounds; they do not
certify a spin-2 thermal seed or a gap.

```sh
julia --startup-file=no --project=models/spin2 models/spin2/test/l2_bounds.jl
julia --startup-file=no --project=models/spin2 models/spin2/scripts/repeated_l2_probe.jl \
  models/spin2/results/repeated-l2-symmetric-replay.json
julia --startup-file=no --project=models/spin2 models/spin2/scripts/loewner_compression_probe.jl \
  models/spin2/results/loewner-compression-replay.json
```

### Faster exact finite-energy verification

The current verifier checks tensor covariance entry by entry using the sparse
ladder generators and multiplies sparse polynomial tensors. Exact covariance
and the checked complete irreducible bases justify highest-weight extraction
by Schur's lemma; optional full-image checks remain available. It also directly
checks commutation with transposition followed by a pi rotation and splits each
transfer block into the two eigenspaces of that involution. Its 92 tests compare
the complete exact norm and energy integers before and after these reductions.
They pass both with one FLINT thread and with four.

Interval contractions now share powers across requested power-of-two lengths
and evaluate traces of products without forming unnecessary dense products.
The measured D=295 matrix-construction stage fell from about 390 seconds to
62 seconds, before the symmetry split and interval contractions. These are
stage timings, not a claim about the complete pipeline's speedup.

For the larger witness, the current command is:

```sh
NEMO_THREADED=1 julia --startup-file=no --project=models/spin2 -e '
  include("models/spin2/scripts/certify_reduced_trial.jl")
  main("models/spin2/results/trial-295-replay.json",
       "models/spin2/results/trial-295-refined.json";
       dimension=295, lengths=(128,256,512,1024), verify_images=false)'
```

Only a result with `complete=true` and strict energy intervals is a completed
finite-energy certificate. Superseded partial runs are not acceptance records.

The completed `results/trial-295-symmetric-certificate.json` gives these strict
finite trial-energy intervals:

| Length | Lower endpoint | Upper endpoint |
| --- | --- | --- |
| 128 | -4.761288829076 | -4.761288829072 |
| 256 | -4.761247250800 | -4.761247250796 |
| 512 | -4.761246333989 | -4.761246333985 |
| 1024 | -4.761246333503 | -4.761246333499 |

Thus **E0(N) < -4.7612463 N** at these four lengths. The largest arithmetic
block is 978, reduced from 1921. The complete run took about 549 seconds with
four FLINT threads. These intervals enclose trial energies; they are not
two-sided ground-energy enclosures and do not certify a gap.

### Local MPO norm and squaring certificates

`CertifiedGaps.CertifiedMPO.Norms` (`julia/src/CertifiedMPO/Norms.jl`) bounds an MPO without expanding the physical Hilbert
space. If the previous partial operators obey `norm(L_b) <= l_b`, and the
next core is `W[b,s,t,a]`, form the nonnegative scalar matrix
`C_a[s,t] = sum_b abs(W[b,s,t,a])*l_b`. The next partial operator norm is at
most `sqrt(norm(C_a,1)*norm(C_a,Inf))`, with induced matrix norms here.
The script uses integer comparisons for all square-root upper bounds and
combines forward and backward witnesses at every cut.

The 33 original checks include the exact norm one of a 128-site identity MPO
and independent dense controls on random MPOs. The ITensors experiment in
`scripts/mpo_norm_probe.jl` concerns **quantized approximate MPOs only**;
it does not certify their distance from the exact thermal operator. The first
four-site beta=1 candidates have bounds above 25 despite numerical
Hilbert--Schmidt norm one. The triangle inequalities lose too much
cancellation. The calculation is scalable, but its present bounds are too
loose for the thermal certificate.

The next construction regards each local core as a map from `(left,input)`
to `(right,output)`. The MPO norm is bounded by the product of these local
operator norms. A numerical diagonal-gauge search reduces this product;
strictly positive rational gauges then cancel exactly on every internal bond.
For each proposed local bound `alpha`, rigorous Arb Cholesky verifies that
`alpha^2 I - W'W` (or the smaller `alpha^2 I - WW'`) is positive definite.
Floating-point singular values only propose bounds; they cannot accept one.

`results/mpo-local-norm-certificate.json` and its saved rational witnesses give:

| Sites | Candidate beta | Initial local norm product | Certified optimized product |
| --- | --- | --- | --- |
| 4 | 1/8 | 3.11769 | 0.16689992 |
| 8 | 1/8 | 38.97081 | 0.04945418 |
| 4 | 1 | 2.94700 | 1.01866711 |
| 8 | 1 | 11.93553 | 1.69278948 |

The candidates have numerical Hilbert--Schmidt norm one. The last two bounds
are still too loose. All four saved certificates pass independent replay
without ITensors or numerical eigensolvers. The norm test file now has 68
passing checks, including intentionally invalid local bounds.

`CertifiedGaps.CertifiedMPO.Squaring` (`julia/src/CertifiedMPO/Squaring.jl`) uses these norms to certify compression of an exact
rational MPO squared. For tail operators `R_a`, define the coefficient map
`R(v)=sum_a v_a R_a`, with Hilbert--Schmidt operator norm at most `gamma`.
Use a separately gauged copy `A_a=g_a R_a`, and let the block row
`[A_1 ... A_D]` have operator norm at most `b`. Then

```
norm(sum_ab V_ab A_a R_b, HS) <= b * gamma * norm(V, HS).
```

This follows by factoring the sum into that block row and the stacked
operators `sum_b V_ab R_b`. It bounds squared tails without a four-layer
contraction. Exact local compression residuals and checked prefix-map norms
then give a telescoping global Hilbert--Schmidt error bound. A Frobenius cap
and optional two-site operator blocks strengthen the tail bounds. The
implementation does not require the approximate MPO to be Hermitian.

At D=8 and initial beta=1/8, the first eight-site relative squaring-error
bound was 6.44554; gauging the first copy reduces it to 1.94875. At the
smaller initial beta=1/128, completed bounds are 0.020357, 0.053945, and
0.134228 for 8, 16, and 32 sites. All three saved rational witnesses pass
independent replay. These are compression errors relative to the computed
output's Hilbert--Schmidt norm, not thermal-trace error bounds. They remain
too large for the required repeated squaring. The initial candidate's
distance from the exact Gibbs amplitude is also still uncertified.

A smaller starting step and separate input/output ranks improve this test.
For the **same exact sixteen-site input**, generated at beta=1/1024 with
input dimension 16, two-site norm blocks and certified compression give:

| Output dimension | Relative squaring-error upper bound |
| --- | --- |
| 16 | 0.00542631 |
| 24 | 0.000693457 |
| 32 | 0.000010764 |
| 48 | 0.000001431 |
| 64 | 0.000000000000481 |

`scripts/refine_mpo_squaring.jl` reuses saved exact input cores and gauges,
so this comparison changes the output rank (and raises coefficient-rounding
precision from 40 to 48 bits for the refinements), without repeating the
thermal evolution. The dimension-64 refinement reaches a relative error
below **4.81e-13**. This remains a **single squaring certificate**, not a
completed repeated-squaring or Gibbs certificate. The completed results are
in `mpo-squaring-small-step-refined.json`, `mpo-squaring-rank-refinement.json`,
and `mpo-squaring-rank-refinement-large.json`. All saved refinements pass
independent replay.

The 59 squaring checks include exact dense controls, gauge invariance,
blocked-core contraction, and rejection of altered residual bounds and
output tensors. Witnesses retain all input and output cores, intermediate
bridges, and local norm bounds. Replay commands are:

```sh
julia --startup-file=no --project=julia julia/test/certified_mpo/mpo_norms.jl
julia --startup-file=no --project=julia julia/test/certified_mpo/mpo_squaring.jl
julia --startup-file=no --project=models/spin2 models/spin2/scripts/verify_mpo_norm.jl \
  models/spin2/results/mpo-local-norm-certificate-n8-step1-witness.json
julia --startup-file=no --project=models/spin2 models/spin2/scripts/verify_mpo_squaring.jl \
  models/spin2/results/mpo-squaring-high-temperature-n32-witness.json
```

```sh
julia --startup-file=no --project=julia julia/test/certified_mpo/positive_bounds.jl
julia --startup-file=no --project=models/spin2 models/spin2/scripts/positive_compression_probe.jl \
  models/spin2/results/positive-compression-weighted-replay.json
```

### A polynomial thermal majorant with a rigorous starting point

`CertifiedGaps.CertifiedMPO.PolynomialBounds` (`julia/src/CertifiedMPO/PolynomialBounds.jl`) supplies a second route to a thermal upper
bound, without a product formula. For `K=H+N*I`, the exact spin-2 bond
spectrum implies `norm(K) <= 5N`. Set `m=2^k`, `delta=beta/(2m)`, and choose
a rational `r >= exp(5N*delta)*(5N*delta)^3/6`. Taylor's remainder bound gives

```
Q = I - delta*K + delta^2*K^2/2 + r*I >= exp(-delta*K) > 0.
Z_N(beta) <= exp(beta*N) * norm(Q^m,HS)^2.
```

Both operators are functions of `K`, so the power comparison is a scalar
eigenvalue argument. Moreover, the exact polynomial trace is at most
`exp(4m*r*exp(5N*delta))` times the true thermal trace. Outward-rounded scalar
checks at **N=256, beta=64** give these upper factors:

| Squarings k | Exact polynomial trace / true trace upper bound |
| --- | --- |
| 24 | 1.177692587 |
| 26 | 1.010236995 |
| 28 | 1.000636180 |
| 30 | 1.000039741 |

Thus 28 squarings cost less than 0.064% in this analytic approximation,
before any MPO compression. The starting object is a specified polynomial
in the exact Hamiltonian. **Its long-chain trace has not yet been contracted
with certified error.** This does not establish the target `Zhat_256(64)<=1.08`.
The 263 tests include the independent exact four-site angular-momentum
spectrum. Results are in `results/polynomial-thermal-targets.json`.

The polynomial representation and physical metric are implemented. For
local basis weights `(1,4,6,4,1)`, a matrix entry's Hilbert--Schmidt weight
is `w_input/w_output`. Local norm checks use rational weighted Gram
inequalities, with outward-rounded Cholesky acceptance. The dyadic norm
endpoints now retain relative accuracy even for exponentially small
normalized operators.

`src/PolynomialMPO.jl` constructs `K` exactly at bond dimension 8
and `Q` at dimension 73. Exact rational row reduction from both ends
reduces the maximum dimension to **34 at N=16 and N=256**, or 33 at N=4,
with **zero initial compression error**. Local intertwining identities
prove equality of the full MPOs; all three portable witnesses pass
independent replay. Results: `results/polynomial-seed-exact.json`.
The earlier rounded dimension-16 seed remains available: its certified
relative error at N=256 is below 1.12e-12, independently replayed.

The exact row-reduced basis needs conditioning before approximate
squaring. `scripts/condition_polynomial_seed.jl` applies invertible
rational matrix changes of basis and saves local equality witnesses.
At N=4, 128-bit gauge proposals and 384-bit ball arithmetic give
right-tail HS bounds below 1.000000001 with zero initial error. Two
dimension-32 squarings then have total certified error **1.225e-13**,
and pass independent replay. The earlier raw exact-basis first step
had an unusable error bound above 7358. Results:
`results/polynomial-seed-conditioned-n4-witness.json` and
`results/polynomial-conditioned-flow-n4.json`.

### Full matrix gauges at 256 sites

The first diagonal-gauge search at N=256 gave a numerical local-norm
product roughly `1.659e75` times the analytic seed bound. A single
Schmidt-based matrix-gauge proposal instead gives about `1.000000675`.
After rational quantization and local interval checks, the rigorous result is

```
norm(B0,operator) <= 1.000000701 * a0
a0 = (1+r+5N*delta+(5N*delta)^2/2)/5^(N/2).
```

This is for the actual saved N=256, dimension-16 starting MPO. The
certificate upper norm is approximately `3.403345322e-90`.
`results/schmidt-seed-norm-n256.json` stores the rational matrix gauges
and local norm witnesses. Independent replay checks the saved source
hash, equality to its physical adjoint, all inverses, and local weighted
Gram inequalities. It does not rerun the Schmidt proposal.

The gauge proposal obtains stable left-prefix Schmidt values by QR/SVD,
then splits their scales between the left and right tensors. Unlike
diagonal scaling, it can rotate the auxiliary basis. Every rational
change `W[j] -> G[j]*W[j]*inv(G[j+1])` cancels exactly in the full MPO.
Candidate generation and rigorous verification remain separate.
Full matrix gauges are supported in the norm, Gram-compression, mixed
product, adjoint-trace, and operator-error residual paths. The last of these
now uses three exact auxiliary-axis contractions to handle noncommuting
matrix gauges, with an interval mode and an exact fallback.

The completed first dimension-16 squaring at N=256 now reaches
`beta=1/2097152` with propagated HS error below **6.762e-8** and an
error-amplification coefficient below **2.001**. The local square and
its verification took about 157 seconds. The earlier diagonal-gauge
run at the same length and dimension had an unusable error bound above
`7.60e25`. These runs also differ in tail grouping and normalization,
so this comparison does not isolate the gauge change alone.
Result: `results/polynomial-schmidt-flow-n256.json`, with its portable
local witness. This is only **one of 28 squarings**: the required beta=64
trace and a uniform spin-2 gap remain unproved.

The first matrix-gauge squaring attempt completed its local checks but
failed during a 192-bit global norm contraction. The flow now obtains
the output HS norm from its verified left-map bounds instead: if those
bounds are `alpha[j]` and the last core is `T`, then
`norm(C,HS) <= norm(T,HS)*prod(alpha)`. A positive rational upper bound
is sufficient for normalization. No dense double-layer contraction is
needed, and the exact recorded scale is retained in the trace formula.
The resumed three-step four-site certificate, with this normalization,
passes independent replay and has total HS error below `1.326e-11`.
Result: `results/polynomial-conditioned-canonical-n4.json`.

Interior bridge contractions can now use outward-rounded intervals to
generate exact rational candidates and residual bounds. Both interval
and exact-rational verifiers accept the small-system fixtures, including
zero and extremely small operators. The final core is still contracted
exactly. The flow saves each raw local witness before normalization;
`resume_from` continues completed steps and `reuse_square_from` recovers
a saved square, while independent replay checks the entire chain.

The optional `error_method=:projection` certifies the cross terms in
the local error sum. For `E=M-U*R`, the witness additionally bounds
`theta >= norm(U'*W*E,HS)`, where `W` is the physical row metric.
Writing `Gamma` for the remaining product-tail norm and `alpha` for
the local prefix norm, the backward suffix-error recurrence is

```
norm(D[j],HS)^2 <= (delta[j]*Gamma[j+1])^2 + alpha[j]^2*e[j+1]^2
                  + 2theta[j]*Gamma[j+1]*e[j+1],   e[N]=0.
```

Here `D[j]` is the actual suffix error. Each new `e[j]` is chosen by outward rounding of the right-hand side's
square root. The verifier checks the projection defect directly; it
does not assume floating-point singular vectors are exactly orthogonal.
`fixed_compression` in the square generator keeps an existing set of
cores and bridges unchanged for a controlled comparison of error bounds.
The flow applies this refinement when `reuse_square_from` is combined
with a different `error_method`. Both rational and interval replay pass
the small-system tests, including a corrupted-bound rejection.

At N=256 and beta=1/2097152, independent replay confirms:

| Retained dimension | Error estimate | Propagated HS error upper bound |
| --- | --- | --- |
| 16 | Triangle | 6.762e-8 |
| 16 | Projection-defect recurrence | 4.259e-9 |
| 32 | Triangle | 8.361e-11 |
| 16 | Projection + rescaled gauges | 7.740e-11 |
| 32 | Projection + 64-bit fitted bridges | 7.413e-12 |

The dimension-16 rows use identical compressed cores and bridges, giving
about a 15.9-fold improvement from the bound alone. Increasing dimension
to 32 instead improves the original bound by about 809-fold. Results:
`results/polynomial-projection-d16-flow-n256.json` and
`results/polynomial-schmidt-d32-flow-n256.json`.
The fourth row retains exactly the original dimension-16 output cores
and normalizing rational, checked in the saved files, and improves its
bound by more than 870-fold. The fifth row uses 64-bit relative rounding
and fits each bridge against the actual rounded left factor; it improves
the original bound by more than 9000-fold. Both pass independent replay.
Results: `results/polynomial-projection-rescaled-d16-n256.json` and
`results/polynomial-projection-d32-fitted-prefix1-n256.json`.
These are first-step error improvements, not a low-temperature trace or
a uniform spin-2 gap certificate.

`rescale_square_gauges` changes `G[j]` to `D[j]*G[j]` with positive rational
diagonal `D[j]`, and divides each bridge's first-input auxiliary index by
the corresponding diagonal entry. The changes cancel exactly, retaining
all compressed cores. Residual and tail bounds must be recertified.
`scripts/gauge_error_screen.jl` numerically screens these rescalings;
its output is explicitly noncertified. `reuse_gauge_power` in the flow
generator requests a full rational/interval certificate for a selected
rescaling of a saved square. Minimizing the initial operator bound alone
need not minimize the propagated compression error.
The accepted rescaling uses Schmidt power 0.2. Its operator-norm bound
is roughly `2.356*a0`, but its smaller compression residual more than
compensates for the increased amplification of the initial error.

The Gram witness now optionally carries a second, independently checked
operator-norm certificate. Residual gauges control the compression tails;
the second gauges can bound the complete input operator more tightly. The
verifier checks both gauge cancellations and all local norm inequalities,
then uses the smaller complete-operator bound in error propagation. It does
not substitute this global bound for any residual tail. The gauge refiner
retains the original operator certificate while choosing new residual gauges.
Small explicit-operator tests cover exact and interval replay, unchanged
compression errors, and rejection of corrupted bounds and boundary gauges.
A four-step quartic integration run passes independent replay, reaching
beta=`1/2048` with HS error at most `1.311e-9`
(`results/quartic-flow-dual-gauge-prefix4-n4.json`). No additional
256-site improvement from this separation has yet been established.

Interval residual replay now retains exact rational gauge matrices and
verifies their inverses exactly, while enclosing the interior transformed
cores in interval arithmetic. This avoids expanding every interior entry
into a large rational coefficient. The terminal core stays exact, and an
inconclusive residual check falls back to exact arithmetic for that core.
The rescaled 256-site certificate passes this independent replay; small
fixtures also pass full exact-rational replay.

New reflected-flow witnesses omit the duplicate `next_candidate` array
and record `next_candidate_encoding="reflected-output"`. The decoder
reconstructs the candidate exactly from the output cores and rational
normalization. Legacy witnesses remain supported. A four-site witness
shrinks from 4,035,669 to 2,800,932 bytes (about 31%); both replay modes and
the adjoint-flow loader accept it and reconstruct the same final candidate.
Result: `results/polynomial-compact-projection-n4.json`.

Exact dyadic coefficients now write large denominators as `2^k`, for example
`7//2^64`. The decoder reconstructs the same integer ratio and continues to
accept every legacy rational token. Reencoding all three witnesses in the
four-site fixture reduces their combined size from 7,562,803 to 5,768,539 bytes
(23.7%), with exact equality of all decoded fields and a successful independent
three-step replay. This reduces storage and parsing input; no end-to-end
speedup is claimed. Result: `results/polynomial-dyadic-projection-n4.json`;
conversion script: `scripts/compact_polynomial_flow.jl`.

Completed witnesses from a partially finished flow can now be resumed. The
generator checks the retained hashes, input chain, errors, and temperatures;
the standalone verifier still replays every retained step. The deliberately
partial four-site fixture resumes to three steps and passes independent replay
(`results/polynomial-partial-resumed-n4.json`). The gauge-error screen now handles
the last step of a multi-step prefix, including its propagated reference norm.

The hybrid 256-site trajectory has completed two independently replayed squarings,
at beta=1/1048576, with propagated error at most `1.747656072e-10`.
The target beta=64 still requires 26 further doublings. Result:
`results/polynomial-hybrid-prefix2-n256.json`.
The trajectory retaining dimension 32 at both steps has also completed
generation, with error at most `1.238239975e-10`; its independent two-step
replay is still running (`results/polynomial-projection-d32-fitted-n256.json`).

### Higher-order initialization through positive squares

`scripts/polynomial_degree_screen.jl` compares exact Taylor majorants. For
`K=H+N*I`, `M=5N`, `delta=beta/(2*2^power)`, and degree `d`, choose

```
r >= exp(delta*M)*(delta*M)^(d+1)/(d+1)!
Q = sum((-delta*K)^k/k!, k=0:d) + r*I.
```

Then `Q >= exp(-delta*K) > 0`, and the exact polynomial partition function
overestimates the Gibbs partition function by at most
`exp(4*2^power*r*exp(delta*M))`. At N=256 and beta=64, certified scalar choices
with less than 0.1% overhead are:

| Polynomial degree | Subsequent squarings | Inflation upper bound, rounded up |
| --- | --- | --- |
| 2 | 28 | 1.000636180 |
| 4 | 21 | 1.000206619 |
| 6 | 18 | 1.000646796 |
| 8 | 17 | 1.000076720 |
| 12 | 16 | 1.000000327 |

These are scalar majorant bounds, not compressed-MPO certificates or thermal
trace estimates. Result: `results/polynomial-degree-screen-n256.json`.

The fourth-order case can reuse the quadratic Gram machinery because

```
P4(x) = (1+x/2+x^2/12)^2 + (x+x^2/2)^2/12 + x^4/72.
```

All coefficients and weights are rational. More generally,
`P_(2k)(x) = integral(exp(-t)*(t+x)^(2k), t=0:infinity)/(2k)!`.
The moment Gram matrix is positive definite; exact rational LDL decomposition
produces a positive sum of squares of degree-k polynomials. The implementation
checks the resulting coefficient identity exactly.

`scripts/quartic_seed_probe.jl` builds three quadratic factors with ITensors,
certifies each Gram compression, and certifies their final weighted sum. If
`sigma=5^(-N/2)`, `A_i=sigma*p_i(-delta*K)`, and
`norm(C_i-A_i'*A_i,HS)<=eta_i`, the assembled seed error is bounded by
`rho + sum(w_i*eta_i)/sigma`, where `rho` certifies final compression.
`scripts/verify_quartic_seed.jl` replays the saved factor witnesses and assembly.
The extra seed construction must be weighed against seven fewer squarings;
no overall runtime improvement or 256-site quartic certificate is claimed yet.
The doubling runner and standalone flow verifier now accept this seed format;
they replay all three factor witnesses before using its error bound. The
four-site seed and a two-step continuation now pass independent replay.
The first four-site prototype (beta=4, power=8, dimension 32) gives a
generator-checked initial error bound of `4.562e-4`. Its final compression
contributes only `1.322e-7`, so the dominant loss is in the factor Gram
certificates. A second small-system run, using the
256-site target's initial step (beta=64, power=21), with factor dimension 34
and output dimension 32, gives a generator-checked bound of `2.072e-8`.
Almost all of that error is in the Gram compression; final assembly contributes
less than `1.759e-15`. Retaining dimension 34 in the products gives
`1.815e-9`, now independently replayed (`quartic-seed-d34-target-n4.json`).

`scripts/refine_quartic_seed.jl` screens positive diagonal rescalings separately
for the three Gram certificates, then recertifies the selected choices. On
that same four-site candidate, powers 0.05, 0.05, and 0.25 give a
independently replayed seed bound of `3.999e-13`. Exact comparisons of the saved
files confirm unchanged inputs, factor outputs, and final candidate, with
more than 4500-fold improvement in the bound. Result:
`results/quartic-seed-d34-refined-n4.json`. The numerical screen itself is
never accepted as proof.

At beta=`1/8192`, two subsequent default-gauge squarings give a replayed
error bound of `3.645e-8` (`results/quartic-flow-prefix2-n4.json`). Applying
the same gauge refinement at each step gives an independently replayed bound of
`2.166e-11`, over 1680 times smaller. Exact file comparisons confirm identical
input cores, output cores, and normalizations at both steps. Result:
`results/quartic-flow-refined-prefix2-n4.json`.
This is still a four-site, high-temperature test; it does not establish a
low-temperature 256-site trace bound.

The optional `exact_inputs=true` seed construction keeps each quadratic factor
exact using two rational compression sweeps and saved local intertwining
identities. The verifier checks those identities before setting its initial
factor error to zero. This may matter at large N: an absolute HS error can
be much larger than a normalized factor's operator norm of order `5^(-N/2)`,
so even small initial rounding errors can be expensive to propagate. This
variant needs conditioning: the unconditioned four-site representation gives
a useless bound above `1.2e4`. Exact right conditioning lowers it to
`1.254e-9`, and subsequent gauge refinement to `3.410e-13`, with zero input
error in all three factors. The latter keeps every factor output and the final
candidate fixed and now passes independent replay. Result:
`results/quartic-seed-conditioned-refined-n4.json`.

The scalar tests include complete four-site spectra, positivity and inflation
checks, and exact sum-of-squares identities through degree 32. The existing
quadratic witnesses keep their original scalar initialization.

### Certified multiplication benchmark

`scripts/benchmark_arb_products.jl` compares the default Arb multiplication with
its certified block routine on synthetic contraction shapes. The largest
products already select the block routine automatically. For shape
`(1024,32)*(32,800)`, the measured times were about 0.50--0.54 seconds versus
0.43--0.44 seconds; the other shapes showed little difference. These are
concurrent-run microbenchmarks, not end-to-end speedups. The verifier's default
multiplication is unchanged. Results: `results/arb-products-benchmark.json`.

```sh
julia --startup-file=no --project=models/spin2 models/spin2/scripts/verify_mpo_norm.jl \
  models/spin2/results/schmidt-seed-norm-n256.json
julia --startup-file=no --project=models/spin2 models/spin2/scripts/verify_polynomial_flow.jl \
  models/spin2/results/polynomial-conditioned-flow-n4.json
julia --startup-file=no --project=models/spin2 models/spin2/scripts/verify_polynomial_flow.jl \
  models/spin2/results/polynomial-schmidt-flow-n256.json
```

`scripts/polynomial_flow_probe.jl` implements Gram products and exact
spatial reflections, allowing alternating compression directions without
an extra rounding step. The completed four-step N=16, D=16 run uses the
earlier rounded seed and reaches only **beta=3.814697265625e-6**. Its
forward normalized HS-error bound is 9.95e-8; backward Gram norm bounds
reduce it to 8.94e-8. Every step passes independent replay. This is a
high-temperature prefix, **not the required beta=64 certificate**.
Results: `results/polynomial-gram-flow-n16.json`.

The exact polynomial's scalar endpoint bound gives a further refinement:
its normalized operator norm is at most
`(1+r+5N*delta+(5N*delta)^2/2)/5^(N/2)`. Propagating this bound through
the recorded positive normalizations reduces the same four-step error
to **4.180e-9**, about 24 times below the original forward bound.
`scripts/polynomial_reference_refinement.jl` independently replays the
source flow before applying this scalar refinement. Result:
`results/polynomial-gram-reference-n16.json`.

`CertifiedGaps.CertifiedMPO.OperatorError` (`julia/src/CertifiedMPO/OperatorError.jl`) separately bounds compression error in
operator norm. Independent gauges on the output prefix and both input
tails enter an exact telescoping sum. For the first saved sixteen-site
step, the bound is **2.673e-18**, compared with the HS-error bound
**8.728e-17**. These control different norms; the new estimate is intended
for the high-entropy part of the flow. Result:
`results/polynomial-gram-operator-error-n16-step1.json`.

The operator-error checker now accepts full matrix gauges on both input
tails and the output prefix. For old first-tail gauges G and new gauges
F, S, and P, the residual transforms as
`E'[q',s,t,a',b'] = sum(P[q',q]*E[q,s,t,a,b]*(G*inv(F))[a,a']*inv(S)[b,b'])`.
The implementation contracts one auxiliary axis at a time, avoiding a dense
Kronecker matrix. It streams the local residuals rather than retaining the
entire chain's residual tensors. Interval replay checks each enclosure and
falls back to an exact local calculation when necessary.

All 138 operator-error tests pass, including exact reconstruction of the
complete physical residual for noncommuting gauges on two-, three-, and
four-site examples, ordinary and Gram products, both replay modes, and
corruption rejection. Unequal auxiliary dimensions are also supported, so
compression of A alone can be certified as the product A times the identity.

The first full-matrix operator-error trial on the second four-site quartic
step gives the independently replayed bound `9.762e-8`, much worse than the
saved HS bound `6.613e-13`. Selecting the prefix and two tail gauge powers
for operator error reduces the former to `2.834e-12`, also independently
replayed. This is a 34,454-fold improvement with the exact input, output,
and original residual bridges fixed, but it is still looser than the HS bound
on this example. The selected powers are `(0.1, 0.05, 0.0)`. Results:
`quartic-operator-error-schmidt-n4-step2.json` and
`quartic-operator-error-refined-n4-step2.json` in `results/`.

`scripts/operator_seed_probe.jl` starts from a verified exact polynomial
MPO, compresses its product with the identity, and saves both HS and operator
error witnesses. `scripts/verify_operator_seed.jl` replays the exact source,
the product compression, and the operator error. At four sites, dimension
16 gives an initial operator bound `1.467e-7`. Gauge refinement reduces it
to `2.280e-13` without changing the candidate or the product witness; both
bounds pass independent replay. Keeping dimension 34 gives an independently replayed
bound `7.289e-17` after gauge refinement, about `1.823e-15` relative to the
analytic operator bound.
Results: `polynomial-operator-seed-refined-n4.json` and
`polynomial-operator-seed-d34-refined-n4.json` in `results/`.
These constructions do not yet provide the missing long-chain operator-error
or trace-moment inputs.

`CertifiedGaps.CertifiedMPO.SchattenBridge` (`julia/src/CertifiedMPO/SchattenBridge.jl`) supplies a conditional transition to a trace
moment: for positive A in dimension d, operator error eta, even p, and
`Tr|B|^p<=M`, the upper bound is
`(1+t)^(p-1)*M + d*(1+1/t)^(p-1)*eta^p` for every positive t.
At d=5^256, p=64, M<=1, and eta<=1e-8 this costs less than 0.05%.
Those long-chain inputs remain to be certified. Its 34 scalar and matrix
checks validate the bridge, not its missing inputs.

### Backward propagation through Schatten norms

The new `SchattenBridge.gram_chain_upper` bounds the desired trace
backward through a compressed Gram chain. Write `D=5^N`, let there be `k`
squarings, and set `p_j=2^(k+1-j)`. If `C_j` approximates
`B_(j-1)'*B_(j-1)` and `B_j` is `C_j/kappa_j` up to spatial reflection,
let `zeta_j` and `sigma_j` be certified HS and operator residual bounds.
Starting from a certified terminal HS bound `m_k`, use

```
m_(j-1) = sqrt_upper(kappa_j*m_j + min(zeta_j, D^(1/p_j)*sigma_j)).
```

An absent operator certificate simply leaves `zeta_j` as the local cost.
The justification is the Schatten triangle inequality and
`norm(B'*B,p)=norm(B,2p)^2`. Therefore `m_0` bounds the initial candidate's
Schatten `2^(k+1)` norm. Add an independently certified initial error,
using the smaller of its HS bound and `D^(1/p_0)` times its operator bound.
The resulting norm bounds the polynomial thermal trace. This recurrence
does not use products of local global-operator norm estimates. At early
steps the dimension factor is close to one; at late steps it can choose
the HS certificate instead.

All square roots, including the dimension factors, are rational upper
bounds checked by exact squaring inequalities. Computing the dimension
factors by repeated square roots avoids constructing integer powers with
hundreds of millions of digits. There are **107 additional passing tests**
on explicit matrix chains, non-Hermitian candidates, spatial reversals,
mixed norm bounds, invalid inputs, and the 28-step exponent at `D=5^256`.
These are tests of the conditional recurrence, not a long-chain trace proof.

`scripts/schatten_flow_refinement.jl` connects the recurrence to saved MPO
flows. It replays the original flow and any supplied local operator-error
witnesses, checks their source hashes, and produces a bound at that flow's
actual inverse temperature. Its `verify_schatten_flow` entry point replays
the saved conclusion independently. A full low-temperature chain remains
necessary for the gap target.

The linear product compressor now also supports interval generation and
replay, with exact fallback for inconclusive residuals and an exact terminal
identity. Its **179 product tests** pass, as do all **138 operator-error tests**.
The rational backend remains available. On the conditioned sixteen-site
seed, the two backends generated exactly the same rounded candidate; this
comparison alone does not certify a speedup or a better error bound.

```sh
julia --startup-file=no --project=julia julia/test/certified_mpo/polynomial_bounds.jl
julia --startup-file=no --project=models/spin2 models/spin2/scripts/polynomial_thermal_targets.jl
julia --startup-file=no --project=models/spin2 models/spin2/test/polynomial_mpo.jl
julia --startup-file=no --project=models/spin2 models/spin2/scripts/verify_polynomial_seed.jl \
  models/spin2/results/polynomial-seed-exact-n256-witness.json
julia --startup-file=no --project=models/spin2 models/spin2/scripts/verify_polynomial_flow.jl \
  models/spin2/results/polynomial-gram-flow-n16.json
julia --startup-file=no --project=julia julia/test/certified_mpo/mpo_operator_error.jl
julia --startup-file=no --project=models/spin2 models/spin2/scripts/verify_mpo_operator_error.jl \
  models/spin2/results/polynomial-gram-operator-error-n16-step1.json
```

### Correcting the trace instead of using only a norm error

For one compression `A -> C` with a certified HS residual `norm(A-C)<=eta`,
the identity

```
norm(A,HS)^2 = 2real(inner(C,A)) - norm(C,HS)^2 + norm(A-C,HS)^2
```

gives an upper certificate with residual contribution `eta^2`. It retains
the signed overlap that is discarded by the triangle bound
`(norm(C,HS)+eta)^2`. For `A=B^2` or `B'*B`, this overlap involves three
MPO layers instead of the four in the direct squared norm.

`scripts/trace_correction_probe.jl` replays each local error certificate,
uses ball arithmetic for the overlap and candidate norm, and checks the
upper bound against a direct uncompressed reference on four sites:

| Output dimension | Triangle upper / reference | Corrected upper / reference |
| --- | --- | --- |
| 16 | 1.155853810 | 1.004971213 |
| 32 | 1.013744177 | 1.000046216 |

At dimension 32 the excess bound shrinks about **297-fold**, with the
same input and compressed candidate. The extra overlap contraction is
essential to that gain. Results: `results/trace-correction-n4-refined.json`.
These are **single-compression certificates for saved candidate MPOs**,
not Gibbs-trace or gap certificates. The multi-step extension below also
accounts for earlier errors.

```sh
julia --startup-file=no --project=models/spin2 models/spin2/scripts/trace_correction_probe.jl \
  models/spin2/results/trace-correction-n4-replay.json \
  models/spin2/results/trace-correction-n4-D16-witness.json \
  models/spin2/results/trace-correction-n4-D32-witness.json
```

### Backward correction across every thermal step

`CertifiedGaps.CertifiedMPO.AdjointTrace` (`julia/src/CertifiedMPO/AdjointTrace.jl`) extends the signed correction to a whole flow.
For reflected Gram steps `A[j+1]=S(A[j]'*A[j])/kappa[j]`, form approximate
backward operators with terminal `G[K]=B[K]` and map
`G[j] ~ B[j]*(S(G[j+1])+S(G[j+1])')/kappa[j]`.
The signed center retains the initial error overlap and every local
compression overlap. Its remaining uncertainty is bounded by

```
2 sum(tau[j]*epsilon[j] + gamma[j]*epsilon[j]^2/kappa[j])
```

where epsilon bounds each primal HS error, tau each backward HS defect,
and gamma the next backward operator norm. The upper endpoint also adds
`epsilon[K]^2`. Exact telescoping proves this formula, including
non-Hermitian candidates. All norm inputs still require certificates.
Backward operator bounds can be propagated as
`g[j] = 2*b[j]*g[j+1]/kappa[j] + tau[j]`, where b bounds the candidate
operator norm, and intersected with separately certified HS bounds.

`scripts/adjoint_trace_probe.jl` implements a complete five-step four-site
benchmark ending at beta=4. Independent replay checks all saved forward
and backward residuals, positive-vector operator bounds, and signed
contractions against an uncompressed Trotter reference:

| Forward final rank | Backward intermediate rank cap | Triangle upper/reference | Corrected upper/reference |
| --- | --- | --- | --- |
| 358 | 625 (uncompressed) | 1.232375999 | 1.018555350 |
| 564 | 625 (uncompressed) | 1.008410682 | 1.000024309 |
| 564 | 358 | 1.008410682 | 1.000301094 |

The second row reduces the excess bound about **346-fold**; compressing
the backward intermediates still gives about **28-fold**. In the last row
the terminal backward operator is the rank-564 forward candidate.
Results: `results/adjoint-trace-n4.json` and
`results/adjoint-trace-compressed-n4.json`. These enclose the complete
four-site **Trotter trace**. Its upper endpoint also bounds the Gibbs trace;
the lower endpoint is not a Gibbs lower bound. No uniform gap follows.

`CertifiedGaps.CertifiedMPO.Products` (`julia/src/CertifiedMPO/Products.jl`) now compresses mixed MPO products of unequal bond
dimensions and contracts three-layer signed overlaps without materializing
the product MPO. Exact local residual witnesses support independent replay.
The streamed first-step N=16 correction reduces the upper bound by
`6.8181e-11` times the candidate squared norm. That comparison is against
the candidate, not an independently computed exact trace.

`scripts/mpo_adjoint_flow.jl` now applies the correction to **all four
steps of the saved N=16 MPO flow**, including its initial seed error.
Both forward and backward MPOs have maximum bond dimension 16. The first
backward attempt without basis conditioning was weaker than the triangle
bound. Reusing the forward rational gauges, reversed for reflection,
removes that loss. Dyadic approximations of the reciprocal gauges keep
the same quality while avoiding expensive denominator growth in bridge
contractions. Every gauge change preserves the represented operator
exactly; rounding the gauge proposal is not an operator approximation.

After independent replay and the backward operator-norm refinement:

| Quantity | Value |
| --- | --- |
| Triangle upper / final candidate squared norm | 1.0000000083597809 |
| Corrected upper / final candidate squared norm | 1.0000000000036466 |
| Corrected upper / corrected lower minus 1 | 6.399e-12 (rounded upward for display) |
| Interval-width improvement using the same corrected lower endpoint | about 1307-fold |

These compare certified enclosures, **not errors against a known exact
trace**. Result: `results/mpo-adjoint-flow-refined-n16.json`. Its source
is `results/mpo-adjoint-flow-dyadic-n16.json`, with four portable product
witnesses. Replay reconstructs every backward input, checks each local
residual and norm bound, and recomputes the three-layer signed overlaps.
The run reaches only **beta=1/262144**, leaving 24 further doublings to
beta=64. The required N=256 low-temperature trace is still missing.

Longer runs now support rigorous upward rounding of the scalar Gram-error
recurrence to 128 relative binary digits. This avoids exponential growth
of exact denominators in otherwise bounded recurrences. A 1000-step
stable recurrence retains denominators at most 150 binary digits, with every
step checked against its exact rational inequality. The new two-step
four-site rounded flow also passes independent replay; old unrounded
certificates remain supported.

The recurrence also uses both independently certified operator norms:
if `norm(A)<=a`, `norm(B)<=b`, and `norm(A-B,HS)<=epsilon`, its linear
coefficient can be the minimum of `2a+epsilon`, `2b+epsilon`, and `a+b`.
The last choice follows directly from
`A'*A-B'*B = A'*(A-B)+(A'-B')*B`. It prevents an artificial
`epsilon^2/kappa` term when normalized operator norms are much smaller
than the available HS error. As a conditional stress test, taking
`a=5^-128`, `b=2a`, `kappa=a`, zero local residual, and `epsilon=1e-12`
gives a next error of `3e-12`; the former coefficient would exceed `1e60`.
This calculation is a propagation stress test; it supplies no local
compression residual. The new formula is implemented in the flow generator and
independent verifier, with a two-step four-site integration certificate
in `results/polynomial-joint-norm-flow-n4.json`.

```sh
julia --startup-file=no --project=julia julia/test/certified_mpo/adjoint_trace.jl
julia --startup-file=no --project=julia julia/test/certified_mpo/mpo_products.jl
julia --startup-file=no --project=models/spin2 models/spin2/test/mpo_adjoint_flow.jl
julia --startup-file=no --project=models/spin2 models/spin2/scripts/adjoint_trace_probe.jl --verify \
  models/spin2/results/adjoint-trace-n4-cap64-witness.json \
  models/spin2/results/adjoint-trace-compressed-n4-cap64-witness.json
julia --startup-file=no --project=models/spin2 models/spin2/scripts/mpo_adjoint_flow.jl --verify \
  models/spin2/results/mpo-adjoint-flow-refined-n16.json
julia --startup-file=no --project=models/spin2 models/spin2/scripts/verify_polynomial_flow.jl \
  models/spin2/results/polynomial-rounded-flow-n4.json
```

### A conditional target requiring one thermal trace

`scripts/one_trace_target.jl` combines temperature lifting and spatial moment
concentration without symmetry-resolved thermal traces. If the shifted trace
`Zhat_n(beta)` is at most `U` and a trial at length `2h*n` proves
`E0(2h*n) < -a*2h*n`, set `D=U^r` and `Q=U^(2h)`. A dominant spatial weight
lower bound `t` yields

\[
S(hn,r\beta)\ge[1+((D-t)/t)^h]^{-2},\qquad
T(2hn,r\beta)\ge[1+(Q-1)^r]^{-2}.
\]

The proof and the exact rule for `t` are in the open
[`optimized-spin1-bounds.tex`](../spin1/experiments/optimized-spin1-bounds.tex).
The 63 tests check the conditional scalar bootstrap and independent rational
moment/Gibbs examples. For `n=128`, `beta=32`, `h=2`, `r=3`, the target `U=1.08`
would give **gap > 0.006039844 for every even L >= 256**. The length-512 trial
condition is certified at `a=4.76122`, but the length-128 trial proves that
the required thermal upper bound is **impossible at this shift**:
`Zhat_128(32) > exp(32*128*(4.761268026385-4.76122)) > 1.1967`.
The exact rational rejection screen is recorded in `one-trace-targets.json`.
This rules out the illustrated seed, not the spin-2 gap. Larger thermal
lengths and more accurate trial shifts must be explored. The scalar gap
values remain explicitly conditional, not measured purities or a gap proof.

With the new D=295 certificate, `results/one-trace-current-targets.json` uses
`n=256`, `beta=64`, `a=4.7612463`, and the certified length-1024 trial. A
rigorous upper trace **Zhat_256(64) <= 1.08** would now give
**gap > 0.003019921 for all even L >= 512**. The known variational trace
lower bound does not reject this target. The thermal inequality itself is
still missing; the JSON keeps `uniform_gap_certified=false`.

## Charge-preserving certificate generation

The opt-in `conserve_charge=true` path now keeps explicit magnetization sectors
through Schmidt gauge balancing, product contractions, and SVD compression.
The shared implementation is `julia/src/CertifiedMPO/ChargeBlocks.jl`
(`CertifiedGaps.CertifiedMPO.ChargeBlocks`); since the move to the core package there is
only the one toolkit, used by every spin-2 script.
The lower-level `product_candidate`, `square_candidate`, and
`schmidt_gauge_proposal` accept `physical_charges=-2:2`. Other spin dimensions
use their own integer charge labels (twice magnetization also works).

Every nonzero must obey `q_right = q_left + m_output - m_input`.
Charge inference checks this exactly, including disconnected support. Dense
gauges can hide this structure without breaking the physical symmetry; such
inputs are rejected by this opt-in path rather than rounded into a different
operator. Each SVD is split into charge blocks, then the largest singular
values across **all** blocks share one `maxdim` budget. Quantization preserves
the structural zeros. The general rational/interval certificate verifiers and
witness formats are unchanged and do not trust the symmetry optimization.

The first saved 16-site benchmark uses the same exactly conditioned seed,
identity second input, no extra gauge, `maxdim=34`, 64 rounding bits,
192-bit interval arithmetic, and one BLAS thread on both paths:

| Measurement | Dense | Charge blocks |
|---|---:|---:|
| Candidate generation | 8.48 s | 5.18 s |
| General certificate replay | 6.89 s | 5.09 s |
| Nonzero output core entries | 261,415 (72.5%) | 32,516 (9.02%) |
| Complete product witness | 11.42 MB | 4.95 MB |
| Certified HS compression error | 3.806e-9 | 2.574e-9 |

These are one full generation and replay per path on a machine with concurrent
work, not a scaling prediction. The error is an upper bound on approximation
error, **not a spectral gap bound**. The error upper bound is not guaranteed to
improve for every input or rank budget.

A separate Gram-product kernel uses a saved physical core and an explicitly
synthetic, charge-compatible bridge. Across three warmed repetitions its median
time falls from 2.962 s to 0.159 s (18.7 times faster), and Julia-tracked allocated
bytes fall from 727 MB to 86 MB. These allocation totals exclude native FLINT
allocations and are not peak resident memory. The full pipeline still stores
dense array shapes and uses the general dense norm and residual checker, so
the kernel improvement is larger than the measured full-run improvement.
Results, input hash, parameters, and witness hashes are saved in
[`charge-block-benchmark-n16.json`](results/charge-block-benchmark-n16.json).

Validation passes 98 new checks through both the legacy and shared entry
points, plus 427 existing product, squaring, and gauge checks. The four-site
seed in `results/polynomial-operator-seed-charge-d34-n4.json` passes its complete
replay. The two-step reflected flow in `results/polynomial-charge-d34-n4.json`
also passes the standalone `verify_polynomial_flow.jl` checker, including the
original exact seed and both saved square witnesses. That integration run
reaches only beta=9.5367431640625e-7, not the target beta=64.

```sh
julia --project=julia julia/test/certified_mpo/charge_blocks.jl
julia --project=models/spin2 models/spin2/scripts/benchmark_charge_blocks.jl \
  models/spin2/results/polynomial-seed-conditioned-n16-witness.json \
  models/spin2/results/charge-block-benchmark-rerun-n16.json
julia --project=models/spin2 -e 'include("models/spin2/scripts/operator_seed_probe.jl"); operator_seed_probe("models/spin2/results/polynomial-seed-conditioned-n4-witness.json", "models/spin2/results/charge-seed-rerun-n4.json"; maxdim=34, conserve_charge=true)'
julia --project=models/spin2 -e 'include("models/spin2/scripts/polynomial_flow_probe.jl"); flow_probe("models/spin2/results/polynomial-seed-conditioned-n4-witness.json", "models/spin2/results/charge-flow-rerun-n4.json"; steps=2, output_dimension=34, gauge_mode=:schmidt, tail_blocksize=1, conserve_charge=true)'
```

The flow currently requires reflected outputs and a fresh compression; the
non-reflected recanonicalizer and reuse of a fixed dense compression are not
enabled in this mode. Existing defaults remain available for replay and
comparison. This is U(1) reduction; full SU(2) multiplet compression and blockwise
interval verification remain further optimization work. No uniform spin-2 gap
is certified by these benchmarks.

## References

- [ITensor custom site types](https://docs.itensor.org/ITensorMPS/stable/examples/Physics.html)
- [ITensor quantum-number DMRG](https://docs.itensor.org/ITensorMPS/stable/tutorials/QN_DMRG.html)
- [ITensor excited-state DMRG](https://docs.itensor.org/ITensorMPS/stable/examples/DMRG.html)
- [Schollwöck and Jolicoeur, spin-2 gap and hidden order](https://arxiv.org/abs/cond-mat/9501115)
  (numerical physical context, not an input certificate)
- [Todo and Kato, Cluster Algorithms for General-S Quantum Spin Systems](https://arxiv.org/abs/cond-mat/9911047), Table I: spin-1 gap 0.41048(6), spin-2 gap 0.08917(4), and spin-2 correlation length 49.49(1). These are numerical benchmarks, not certificate inputs.
