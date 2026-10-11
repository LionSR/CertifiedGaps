# Seed search from exact diagonalisation for two models of the survey:
#   * spin-½ Heisenberg ring with staggered field h_s = 1 (models/heisenberg_staggered_field):
#     rings up to L = 16, sectors above `dense_max` states by Lanczos;
#   * the spin-1 Heisenberg ring of the note (models/spin1_haldane_dmrg): rings up to L = 12,
#     large sectors by Lanczos, where the survey found no seed with n₀ ≤ 6.
# Run: julia --project=julia julia/scripts/ed_seeds.jl

using ThermalMPOGap, Printf, LinearAlgebra
using ThermalMPOGap.ExactDiagonalization, ThermalMPOGap.Models, ThermalMPOGap.Seeds

BLAS.set_num_threads(1)   # one BLAS thread per process; oversubscription is far slower

function survey(name, model, lengths, ns, βs; kwargs...)
    println("\n## $name")
    spectra = Dict{Int,PartialSpectrum}()
    for L in lengths
        t = @elapsed spectra[L] = spectrum(model.chain, L; kwargs...)
        s = spectra[L]
        @printf("L = %2d: dim %8d, exact %-5s, E₀ = %.10f, gap %.6f  (%.1f s)\n",
            L, s.N, s.exact, s.E[1], s.E[2] - s.E[1], t)
        flush(stdout)
    end
    results = scan_seeds(spectral_defects(spectra), ns, βs)
    for n in ns
        best = maximum(r -> min(1 - r.εₛ, 1 - r.εₜ), filter(r -> r.n₀ == n, results))
        @printf("n₀ = %d: max over β₀ of min(Pₛ, P_τ) ≥ %.6f\n", n, best)
    end
    best = best_seed(results)
    best === nothing ? println("criterion not met on this grid") :
        @printf("best seed (n₀, β₀) = (%d, %g): εₛ = %.4e, εₜ = %.4e, Δ_L ≥ %.6f\n",
            best.n₀, best.β₀, best.εₛ, best.εₜ, best.bound)
end

βgrid = 0.25:0.25:32
survey("Staggered-field Heisenberg, h_s = 1", staggered_heisenberg(1.0),
    [2, 4, 6, 8, 12, 16], [2, 4, 6, 8], βgrid; dense_max = 3000, k = 80)
survey("Spin-1 Heisenberg ring (the note's model)", heisenberg(1),
    [2, 4, 6, 8, 10, 12], [2, 4, 6], βgrid; dense_max = 3000, k = 80)
