# Checker for the cluster-chain certificates.
#
#     julia --project=code/julia code/models/cluster_chain/check.jl          # exact parts
#     julia --project=code/julia code/models/cluster_chain/check.jl --full   # + log Z
#
# For every point of POINTS (src/common.jl), `check_point` requires the model-local source
# hashes, the point, seed and exact RP record to be the current definitions, and reruns the
# interval bound and the exact chain from the stored εₛ, εₜ. `--full` also recomputes the
# closed-form log Z enclosures and requires the fresh defects to be at most the stored ones.
# Finally results/summary.json must list the certificates' rows.
include(joinpath(@__DIR__, "src", "common.jl"))
full = "--full" in ARGS
for p in POINTS
    data = read_result(certificate_path(p.g), Dict{String,Any})
    seconds = @elapsed b, c = check_point(data, p.g, p.n0, p.beta0; full)
    println("g = $(p.g): Δ ≥ $b (chain > $(Float64(c))) for even L ≥ $(p.n0) cells; checked",
        full ? " with log Z" : "", " in $(round(seconds; digits = 1)) s")
end
summary = read_result(joinpath(ROOT, "results", "summary.json"), Dict{String,Any})
[(r["g"], r["n0_cells"], r["eps_s"], r["eps_t"]) for r in summary["points"]] ==
    [let d = read_result(certificate_path(p.g), Dict{String,Any})
         (d["g"], d["seed"]["n0_cells"], d["seed"]["eps_s"], d["seed"]["eps_t"])
     end for p in POINTS] || error("results/summary.json does not match the certificates")
println("all $(length(POINTS)) certificates and the summary verified")
