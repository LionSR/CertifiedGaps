# Routes considered

The target is an exact RP identity and a verified full-spectrum seed at strong rung coupling.

* Dense ED with magnetization sectors suffices for 4 rungs, but the largest 8-rung
  sector has dimension 12,870 and the 10-rung sector 184,756. Use translation by
  one rung, leg exchange, and spin inversion to reduce numerical scan blocks.
  Dense cost is cubic in block size and quadratic in storage; beta scans are cheap
  once all eigenvalues have been computed. No Lanczos truncation is a certificate.
* For the rigorous 4/8-rung seed, residual and orthogonality norm bounds on a full
  dense eigendecomposition give simultaneous eigenvalue enclosures. Rigorous
  floating-point dot-product error bounds can avoid cubic BigFloat arithmetic.
* Strong-rung expansions could reach much larger rings, but a useful explicit
  convergent remainder is substantially more work than the requested ED seed.
* A product-rung perturbation bound is simple but too loose for these purities.
* Transfer-operator contraction is attractive beyond ED but requires an additional
  analytic estimate not supplied by the core.

Chosen route: symmetry-resolved full dense numerical scan, then verified full
spectra for the best admissible even seed at J_perp=4. Only even n0 are valid in
Appendices D–F; odd-n scans, if made, are diagnostics, not certificate candidates.
The threshold means success on the stated finite grid, not a phase transition or
proof of absence of seeds at other sizes/temperatures.

Environment: initial approved juliaup installer request returned HTTP 403.

Resource observation: running six numerical parameters concurrently with dense
complex intermediate matrices exceeded the 17-GiB memory limit during a
10-rung block. The J=7/2 process was killed (exit 137; cgroup oom_kill=1),
while all completed small spectra and the J=4 certificate persisted. A tested
sparse-to-real change of basis (`sparse_scan_basis.jl`) avoids those large
complex dense intermediates; use it for retries and future numerical scans.
Do not infer missing large-ring results from an interrupted process.
