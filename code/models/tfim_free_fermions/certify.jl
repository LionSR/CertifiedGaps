# Generate the free-fermion TFIM certificates (one per point of POINTS in src/common.jl).
#
#     julia --project=code/julia code/models/tfim_free_fermions/certify.jl
#
# Each certificate holds the interval enclosures of log Z_{n₀}(β₀), log Z_{2n₀}(β₀) and
# log Z_{2n₀}(2β₀) from the closed form, the rational seed defects, the interval bound
# `certified_gap_bound` and the companion exact chain `rational_gap_chain`.
include(joinpath(@__DIR__, "src", "common.jl"))

rows = Any[]
for p in POINTS
    seconds = @elapsed data = certify_point(p.g, p.n0, p.beta0)
    write_result(certificate_path(p.g), data)
    s = data["seed"]
    push!(rows, Dict{String,Any}("g" => data["g"], "h" => "0", "n0" => p.n0,
        "beta0" => s["beta0"], "eps_s" => s["eps_s"], "eps_t" => s["eps_t"],
        "interval_bound_decimal" => s["interval_bound_decimal"],
        "chain_bound_decimal" => data["chain"]["gap_strict_lower_decimal"],
        "exact_gap" => data["exact_gap"]))
    println("g = $(p.g): Δ ≥ $(s["interval_bound_decimal"]) (chain ",
        data["chain"]["gap_strict_lower_decimal"], ") for even L ≥ $(p.n0), ",
        "in $(round(seconds; digits = 1)) s")
end
write_result(joinpath(ROOT, "results", "summary.json"),
    Dict{String,Any}("format" => FORMAT * " summary", "points" => rows))
