# Checker for the staggered-field Heisenberg certificates (core format, one per field).
#
#     julia --project=julia models/heisenberg_staggered_field/check.jl          # exact parts
#     julia --project=julia models/heisenberg_staggered_field/check.jl --full   # + spectra
#
# For every point, `check_point` requires the model-local source hashes, the stored model, RP
# data, seed, blocking and params to be the point's definitions in src/common.jl (not merely
# self-consistent), the stored εₛ, εₜ to be at most the seed defects the Lean theorem assumes
# (`lean_points`), and then reruns the exact chain and the interval bound. `--full` also
# rebuilds and re-verifies the 8- and 16-site spectra of every point and requires the fresh
# defects to be at most the stored εₛ, εₜ. Finally results/summary.json must list exactly the
# certificates' rows.
include(joinpath(@__DIR__, "src", "common.jl"))
full = "--full" in ARGS
lean = Dict(p.params => p for p in lean_points(LEAN))
for (h, β0) in SEEDS
    seconds = @elapsed check_point(certificate_path(h), site_model(h); point(h)...,
        lean = lean[[h, β0]], spectra = full)
    println("h = $h: certificate checked", full ? " with spectra" : "",
        " in $(round(seconds; digits = 1)) s")
end
summary = read_result(joinpath(ROOT, "results", "summary.json"), Dict{String,Any})
summary["points"] ==
    [certificate_row(read_result(certificate_path(h), Dict{String,Any})) for h in FIELDS] ||
    error("results/summary.json does not match the certificates")
println("all $(length(SEEDS)) certificates and the summary verified")
