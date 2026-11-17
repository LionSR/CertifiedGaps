# Checker for the XY-in-a-field certificates.
#
#     julia --project=code/julia code/models/xy_field/check.jl          # exact parts
#     julia --project=code/julia code/models/xy_field/check.jl --full   # + log Z
#
# For every point of POINTS (src/common.jl), `check_point` requires the model-local source
# hashes, the parameters, seed and exact RP record to be the current definitions, and reruns the
# interval bound and the exact chain from the stored εₛ, εₜ. `--full` also recomputes the
# closed-form log Z enclosures and requires the fresh defects to be at most the stored ones.
# Finally results/summary.json must list the certificates' rows.
include(joinpath(@__DIR__, "src", "common.jl"))
full = "--full" in ARGS
for p in POINTS
    data = read_result(certificate_path(p), Dict{String,Any})
    seconds = @elapsed b, c = check_point(data, p; full)
    println("(Jx, Jy, h) = ($(p.Jx), $(p.Jy), $(p.h)): Δ ≥ $b (chain > $(Float64(c))) ",
        "for even L ≥ $(p.n0); checked", full ? " with log Z" : "",
        " in $(round(seconds; digits = 1)) s")
end
summary = read_result(joinpath(ROOT, "results", "summary.json"), Dict{String,Any})
[(r["Jx"], r["Jy"], r["h"], r["n0"], r["eps_s"], r["eps_t"]) for r in summary["points"]] ==
    [let d = read_result(certificate_path(p), Dict{String,Any})
         (d["Jx"], d["Jy"], d["h"], d["seed"]["n0"], d["seed"]["eps_s"], d["seed"]["eps_t"])
     end for p in POINTS] || error("results/summary.json does not match the certificates")
println("all $(length(POINTS)) certificates and the summary verified")
