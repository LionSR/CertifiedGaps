using JSON, Printf
rows=JSON.parsefile(joinpath(@__DIR__,"costs.json"))["rows"]
open(joinpath(@__DIR__,"TABLES.md"),"w") do io
    println(io,"# All 43 seed targets\n\nNON-CERTIFYING spectral counts and runtime projections. Dimensions and target constants are source-derived.\n")
    println(io,"`η = min(10⁻⁶, min(εₛ−εₛ,source, εₜ−εₜ,source)/16)`; cutoff `E₀+(log D+log(1/η))/β`. `r=η/(4β)` is a conservative planning radius. Gram minutes exclude vector generation, anchors, analytic proofs, and encoding costs.\n")
    println(io,"| Model / parameter | D | β | εₛ | εₜ | η | r | cutoff − E₀ | k | Gram minutes, balanced–concentrated |")
    println(io,"|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|")
    for r in rows
        @printf(io,"| %s | %d | %.6g | %.9g | %.9g | %.3g | %.3g | %.5g | %d | %.3g–%.3g |\n",
            replace(r["model"],"//"=>"/"),r["dimension"],r["beta"],r["eps_s"],r["eps_t"],
            r["tail_budget"],r["uniform_eigen_radius_budget"],r["cutoff_above_ground"],
            r["levels_below_cutoff"],60r["gram_hours_balanced"],60r["gram_hours_concentrated"])
    end
end
