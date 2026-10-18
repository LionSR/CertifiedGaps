# Rational data for `CertifiedGaps.OpenCertificate` (lean/CertifiedGaps/Open/Certificate.lean).
#
# Given the seed half-lengths N₀ (ring: P_s(2N₀, β₀), P_τ(4N₀, β₀) ≥ 1 - u₀) and K₀ (boundary:
# e_{2K₀}(β₀) ≤ E₀), the open-chain seed P_τ(4N₀, β₀) ≥ 1 - η₀ and β₀ (exact rationals), this
# computes short decimal upper bounds along the base steps of the boundary induction,
#     εᵢ = (C u₀)^{2^i} / C,                           C ≥ G(u₀) = ½((1-u)⁻¹ + (1-u)⁻² + (1-u)⁻³)²,
#     qᵢ^{bᵢ/g} ≥ d(εᵢ₊₁)^{aᵢ/g},                       d(ε) = ε / (2 - 3ε),
#     Eᵢ₊₁ ≥ F(qᵢ, Eᵢ, 1 - ηᵢ),                         F(q, e, t) = q/(1-q) (1/((1-e)² t) - 1),
#     ηᵢ₊₁ ≥ f(1 - (1 - ηᵢ)(1 - εᵢ)(1 - εᵢ₊₁)(1 - Eᵢ₊₁)),   f(z) = z² / (2(1 - z)²),
# with aᵢ = 4Nᵢ - 2Kᵢ, bᵢ = 4Nᵢ, g = gcd(aᵢ, bᵢ), Nᵢ = N₀ 2^i, K₀ given and Kᵢ = Nᵢ for i ≥ 1.
# For each number of steps m it hands over at x = max(Eₘ, ηₘ, εₘ) ≤ 1/100, which gives
#     Δ_n > log(1/(10x)) / (β₀ 2^m)   for every odd n ≥ 4 N₀ 2^m,
# and prints every m with its bound and the certificate with the best bound. Lean re-checks every
# inequality by `decide +kernel`; this script only searches, in exact `Rational{BigInt}`
# arithmetic (the root d^{a/b} is taken in `BigFloat` and then checked exactly).
#
# Usage: julia lean/scripts/open_certificate.jl N₀ K₀ u₀ E₀ η₀ β₀ [mmax] [digits]
#   e.g. julia lean/scripts/open_certificate.jl 60 71 0.052 0.14 0.21 49
# (the seeds of the note's eq. (C1)). Arguments are parsed as literals, never evaluated.

include(joinpath(@__DIR__, "chain_certificate.jl"))

coefficient(u) = ((1 - u)^-1 + (1 - u)^-2 + (1 - u)^-3)^2 / 2
tail(ε) = ε / (2 - 3ε)
overlap(q, e, t) = q / (1 - q) * (1 / ((1 - e)^2 * t) - 1)
boundary_half(N₀, K₀, i) = i == 0 ? K₀ : N₀ * 2^i

"Smallest `digits`-digit decimal `q` with `q^b ≥ d^a`, found from the real root `d^(a/b)`."
function suppression(d, a, b, digits)
    a == 0 && return one(d)
    q = round_up(Rational{BigInt}(BigFloat(d)^(BigFloat(a) / b)), digits)
    k = digits - 1 - floor(Int, log10(BigFloat(q)))
    ulp = big(1) // big(10)^max(k, 0) * big(10)^max(-k, 0)
    while q^b < d^a
        q += ulp
    end
    return q
end

"The base steps up to `mmax`, stopping early if a defect leaves the admissible range."
function steps(N₀, K₀, u₀, C, E₀, η₀, mmax, digits)
    ε(i) = (C * u₀)^(big(2)^i) / C
    q, E, η = Rational{BigInt}[], [E₀], [η₀]
    for i in 0:(mmax - 1)
        len = 4 * N₀ * 2^i
        a = len - 2 * boundary_half(N₀, K₀, i)
        g = gcd(a, len)
        qi = suppression(tail(ε(i + 1)), a ÷ g, len ÷ g, digits)
        (qi < 1 && E[end] < 1 && η[end] < 1) || break
        Ei = round_up(overlap(qi, E[end], 1 - η[end]), digits)
        Ei <= 1 || break
        prod = (1 - η[end]) * (1 - ε(i)) * (1 - ε(i + 1)) * (1 - Ei)
        prod > 0 || break
        push!(q, qi)
        push!(E, Ei)
        push!(η, round_up(f(1 - prod), digits))
    end
    return q, E, η, ε
end

function open_main(args)
    N₀, K₀ = parse.(Int, args[1:2])
    u₀, E₀, η₀, β₀ = parse_rational.(args[3:6])
    mmax = length(args) ≥ 7 ? parse(Int, args[7]) : 8
    digits = length(args) ≥ 8 ? parse(Int, args[8]) : 4
    C = round_up(coefficient(u₀), digits)
    (C <= 6 && C * u₀ <= 1 && u₀ <= 1 // 10) ||
        error("ring seed u₀ = $(Float64(u₀)) is too large: G(u₀) ≈ $(Float64(C))")
    q, E, η, ε = steps(N₀, K₀, u₀, C, round_up(E₀, digits), round_up(η₀, digits), mmax, digits)
    best = nothing
    println("-- m   x            gap bound        odd n ≥")
    for m in 1:length(q)
        x = round_up(max(E[m + 1], η[m + 1], ε(m)), digits)
        x <= 1 // 100 || continue
        bound = log(1 / (10 * BigFloat(x))) / (BigFloat(β₀) * 2^m)
        println("-- ", rpad(m, 4), rpad(Float64(x), 13), rpad(Float64(bound), 17), 4 * N₀ * 2^m)
        if best === nothing || bound > best[3]
            best = (m, x, bound)
        end
    end
    best === nothing && error("no hand-over x ≤ 1/100 within $mmax steps")
    m, x, bound = best
    println("-- best: m = $m, gap bound log(1/(10x)) / (β₀ 2^m) = ", Float64(bound))
    println("  N₀ := ", N₀)
    println("  K₀ := ", K₀)
    println("  u₀ := ", lean(u₀))
    println("  C := ", lean(C))
    println(lean_list("q", q[1:m]))
    println(lean_list("E", E[1:(m + 1)]))
    println(lean_list("η", η[1:(m + 1)]))
    println("  x := ", lean(x))
end

if abspath(PROGRAM_FILE) == (@__FILE__)
    open_main(ARGS)
end
