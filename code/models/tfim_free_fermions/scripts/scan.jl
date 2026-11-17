# Seed scan (a diagnostic, not a certificate): for each g, the bound
# certified_gap_bound(εₛ, εₜ, β₀) over even n₀ and a grid of β₀, with the defects from the
# interval closed form of src/FreeFermions.jl at reduced precision. Float64 cannot be used:
# log Z ≈ 10³ while the defects are far below its rounding error.
#
#     julia --project=code/julia code/models/tfim_free_fermions/scripts/scan.jl g ns βs
#
# e.g. `scan.jl 2 16:16:256 4:1:64`; g, β₀ are read as rationals (`5//4`).
using CertifiedGaps, Logging
# Failed seeds make the interval recurrence emit NaI warnings; the scan only reads `bound`.
Logging.disable_logging(Logging.Warn)
include(joinpath(@__DIR__, "..", "src", "FreeFermions.jl"))

parse_q(s) = Rational{BigInt}(eval(Meta.parse(s)))
g = parse_q(ARGS[1])
ns = eval(Meta.parse(ARGS[2]))
βs = Rational{BigInt}.(eval(Meta.parse(ARGS[3])))
println("g = $g (exact gap $(2(g - 1)))")
best = nothing
for n in ns
    rowbest = nothing
    for β in βs
        l = (logZ(n, g, β; precision = 256), logZ(2n, g, β; precision = 256),
            logZ(2n, g, 2β; precision = 256))
        d = CertifiedGaps.Seeds.seed_defects(l...)
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
