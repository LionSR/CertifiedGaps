# Example: transverse-field Ising ring seeds from the free-fermion formula
# (`Models.tfim_seed_defects`, survey Proposition 4.13), at g = 2 and g = 1.25, the TFIM open
# targets listed in AGENTS.md (no model directory yet). The BigFloat defects are padded
# heuristically, so the table is screening, not a certificate.
# Run: julia --project=julia julia/examples/tfim.jl

using CertifiedGaps, Printf
using CertifiedGaps.Models: tfim_seed_defects

# BigFloat defects at 2048 bits, padded upward by a relative 1e-30 before rounding to Float64.
padded(ε) = Float64(ε * (1 + big(10)^-30), RoundUp)

function report(g, seeds)
    @printf("\n## g = %s (exact gap 2(g − 1) = %s)\n", g, 2(g - 1))
    println("| n₀ | β₀ | εₛ | εₜ | bound |")
    println("|---|---|---|---|---|")
    for (n, β) in seeds
        εₛ, εₜ = padded.(tfim_seed_defects(n, β, g))
        b = certified_gap_bound(εₛ, εₜ, β).bound
        @printf("| %d | %g | %.4e | %.4e | %s |\n", n, β, εₛ, εₜ,
            b === nothing ? "none" : @sprintf("%.6f", b))
    end
end

report(2.0, [(8, 3), (8, 64), (16, 4), (32, 8), (64, 16), (128, 32), (256, 64)])
report(1.25, [(144, 64), (256, 64)])
