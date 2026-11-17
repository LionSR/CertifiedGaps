# Seed scan (a diagnostic, not a certificate): certified_gap_bound(εₛ, εₜ, β₀) over even n₀ and
# a grid of β₀, with the defects from the interval closed form of src/FreeFermions.jl at 256 bits
# (Float64 cannot resolve defects far below its rounding error of log Z).
#
#     julia --project=code/julia code/models/xy_field/scripts/scan.jl Jx Jy h ns βs
#
# e.g. `scan.jl -1 1//2 1 "[16,32,64]" "2:1//2:40"`; parameters and β₀ are read as rationals.
using CertifiedGaps, Logging
# Failed seeds make the interval recurrence emit NaI warnings; the scan only reads `bound`.
Logging.disable_logging(Logging.Warn)
include(joinpath(@__DIR__, "..", "src", "FreeFermions.jl"))

q(s) = Rational{BigInt}(eval(Meta.parse(s)))
Jx, Jy, h = q(ARGS[1]), q(ARGS[2]), q(ARGS[3])
ns = eval(Meta.parse(ARGS[4]))
βs = Rational{BigInt}.(eval(Meta.parse(ARGS[5])))
println("(J_x, J_y, h) = ($Jx, $Jy, $h), exact gap $(2(h - abs(Jx + Jy)))")
best = nothing
for n in ns
    rowbest = nothing
    for β in βs
        l = setprecision(BigFloat, 256) do
            (logZ(n, Jx, Jy, h, β; precision = 256), logZ(2n, Jx, Jy, h, β; precision = 256),
                logZ(2n, Jx, Jy, h, 2β; precision = 256))
        end
        d = setprecision(() -> CertifiedGaps.Seeds.seed_defects(l...), BigFloat, 256)
        (0 <= d.eps_s < 1 && 0 <= d.eps_t < 1) || continue
        b = certified_gap_bound(d.eps_s, d.eps_t, β).bound
        b === nothing && continue
        row = (n, β, Float64(d.eps_s), Float64(d.eps_t), b)
        (rowbest === nothing || b > rowbest[end]) && (rowbest = row)
    end
    println("  ", rowbest)
    rowbest !== nothing && (best === nothing || rowbest[end] > best[end]) && (global best = rowbest)
end
println("best: ", best)
