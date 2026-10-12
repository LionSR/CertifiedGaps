# Constants of the doubling recurrence (survey §2.2, §2.5, §2.6), recomputed independently.
# Run: julia --project=julia julia/scripts/constants.jl

using ThermalMPOGap, IntervalArithmetic, Printf

B(x) = interval(BigFloat, x)
lower(x) = Float64(inf(x), RoundDown)

println("## Diagonal threshold (Theorem 2.17)")
@printf("ε⋆ = %.12f; exact rational bracket 0.12731 < ε⋆ < 0.12732 holds: %s\n",
    Float64(diagonal_fixed_point()), fixed_point_bracket_holds())

println("\n## The note's seeds εₛ = εₜ = 1/100 (certified, outward rounded)")
c = certified_gap_bound(1 // 100, 1 // 100, 1)
@printf("β₀ Δ_L ≥ %.10f   (note: log 20 = %.10f)\n", c.bound, log(20))
@printf("spin-1 ring, β₀ = 784: Δ_L ≥ %.7f for even L ≥ 2304   (note: %.7f)\n",
    c.bound / 784, log(20) / 784)

println("\n## Odd lengths (Proposition 2.13): β₀ Δ_L ≥ log((1 − 2u₀)/(2u₀)), u₀ = 1 − 0.99³")
u₀ = 1 - (B(99) / 100)^3
odd = log((1 - 2u₀) / (2u₀))
@printf("β₀ Δ_L ≥ %.7f; spin-1: Δ_L ≥ %.7f for odd L ≥ 2305\n", lower(odd), lower(odd / 784))

println("\n## Closed forms (Theorem 2.22) from the majorant recursion, J = 40")
setprecision(BigFloat, 256) do
    for (label, R₀, P₀) in (("(i)  diagonal ε₀ ≤ 1/100", B(3) / 97, B(1) / 3),
                            ("(ii) r₀ ≤ 0.03", B(3) / 100, B(1)),
                            ("(ii) a₀ ≤ 1/100", B(1) / 99, B(1)))
        Λ, ok = majorant_loss(R₀, P₀; J = 40)
        @printf("%-26s Λ₄₀ ≤ %.6f, e^Λ₄₀ ≤ %.4f, tail hypothesis holds: %s\n",
            label, Float64(sup(Λ), RoundUp), Float64(sup(exp(Λ)), RoundUp), ok)
    end
    Λ, _ = majorant_loss(B(3) / 97, B(1) / 3; J = 40)
    @printf("(i) at ε₀ = 1/100: β₀ Δ_L ≥ log(0.97/(3ε₀)) − Λ₄₀ = %.4f, i.e. log(%.5f/ε₀)\n",
        lower(log(B(97) / 3) - Λ), lower(exp(log(B(97) / 300) - Λ)))
end
