# Spin-1 bilinear–biquadratic chain: route analysis

Each unit cell is a single spin-1 site ($d = 3$). The single-cell reflection theorem requires an even seed length $n_0$; odd $n_0$ values serve only as diagnostics.

## Candidate routes

1. **Closed forms.** No closed-form partition function exists for the interacting bilinear–biquadratic chain in the Haldane regime.
2. **Convergent cluster expansions.** Near the isotropic and negative-$\theta$ Haldane regime, there is no trivial product-state limit around which to expand without introducing an auxiliary anisotropy or dimerization parameter.
3. **Translation, magnetization, and $\mathrm{SU}(2)$ symmetry reduction (chosen for finite-size survey).** Combining translation momentum, magnetization, spin flip, and spatial reflection block-diagonalizes complete rings through $L = 12$. To reach $L = 14$ ($3^{14} = 4{,}782{,}969$ states) within 16 GiB of RAM, we further project each $(M = S, k)$ block onto its $\mathrm{SU}(2)$ highest-weight subspace (kernel of $S^+$). Specifically, in sector $(M = S, k)$, the operator $D^\dagger D = S^- S^+$ has exact eigenvalues $J(J+1) - S(S+1)$ for $J = S, \dots, L$. Applying the polynomial projector $\prod_{J > S} \bigl(I - D^\dagger D / [J(J+1) - S(S+1)]\bigr)$ to deterministic Gaussian columns and orthonormalizing via `CholeskyQR2` yields an orthonormal highest-weight basis of dimension $r = d(M = S, k) - d(M = S+1, k)$, which is then diagonalized densely with multiplet weight $2S + 1$.

## Outcome

**No certificate at $n_0 \le 6$.** While the exact reflection-positivity identity holds across $-\pi/4 < \theta \le 0$, complete spectra up to $L = 14$ show that no even seed $n_0 \le 6$ approaches the admissible region at any of the six tested angles (the best purity is $\min(P_{\mathrm{s}}, P_\tau) = 0.329$ at $\theta = 0, n_0 = 6$, compared to the $\approx 0.873$ threshold). Reaching $n_0 = 8$ ($L = 16$) would require diagonalizing a $59{,}143 \times 59{,}143$ $\mathrm{SU}(2)$ block (26.1 GiB per dense matrix), and finite-size scaling of the ground-energy mismatch $\delta(n) = E_0(2n) - 2E_0(n)$ estimates that admissible seeds require $n_0 \approx 20\text{–}48$. Certifying the negative-$\theta$ BLBQ chain therefore requires a large-$n_0$ non-ED method such as polynomial-filtered trace bounds combined with MPS variational upper bounds (as in [`models/spin1/`](../spin1/README.md)).
