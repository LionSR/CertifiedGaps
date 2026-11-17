# Seed scan (a diagnostic, not a certificate): certified_gap_bound(εₛ, εₜ, β₀) over even n₀ and
# a grid of β₀, with the defects from the interval closed form of src/FreeFermions.jl at 256 bits
# (Float64 cannot resolve defects far below its rounding error of log Z).
#
#     julia --project=code/julia code/models/cluster_chain/scripts/scan.jl g ns βs
#
# e.g. `scan.jl 1//2 "[32,64,128]" "8:1:80"`; n₀ counts two-site cells, g and β₀ are rationals.
using CertifiedGaps, Logging
# Failed seeds make the interval recurrence emit NaI warnings; the scan only reads `bound`.
Logging.disable_logging(Logging.Warn)
include(joinpath(@__DIR__, "..", "src", "FreeFermions.jl"))

q(s) = Rational{BigInt}(eval(Meta.parse(s)))
g = q(ARGS[1])
ns = eval(Meta.parse(ARGS[2]))
βs = Rational{BigInt}.(eval(Meta.parse(ARGS[3])))
println("g = $g, exact gap $(2abs(1 - g))")
best = nothing
for n in ns
    rowbest = nothing
    for β in βs
        l = setprecision(BigFloat, 256) do
            (logZ(n, g, β; precision = 256), logZ(2n, g, β; precision = 256),
                logZ(2n, g, 2β; precision = 256))
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
