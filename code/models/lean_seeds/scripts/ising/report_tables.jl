# Print the 43-row NON-CERTIFYING cost table (results/ising/costs.json) as Markdown (METHOD.md).
using JSON, Printf
rows=JSON.parsefile(joinpath(@__DIR__,"..","..","results","ising","costs.json"))["rows"]
println("| Model / parameter | D | β | εₛ | εₜ | η | r | cutoff − E₀ | k | Gram minutes, balanced–concentrated |")
println("|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|")
for r in rows
    @printf("| %s | %d | %.6g | %.9g | %.9g | %.3g | %.3g | %.5g | %d | %.3g–%.3g |\n",
        replace(r["model"],"//"=>"/"),r["dimension"],r["beta"],r["eps_s"],r["eps_t"],
        r["tail_budget"],r["uniform_eigen_radius_budget"],r["cutoff_above_ground"],
        r["levels_below_cutoff"],60r["gram_hours_balanced"],60r["gram_hours_concentrated"])
end
