# Certifies the staggered-field points with the shared core: for each h in `SEEDS`
# (src/common.jl), `GapCertificates.certify_gap` on the one-site rotated-frame model at
# (n₀, β₀) = (8, β₀), with the seed spectra split by the two-site U(1) cell charge `BLOCKING`.
# Writes results/certificate_<h>.json and results/summary.json.
#     julia --project=julia models/heisenberg_staggered_field/certify.jl [h]   (default: all five)
include(joinpath(@__DIR__, "src", "common.jl"))
fields = isempty(ARGS) ? FIELDS : [parse_rational(ARGS[1])]
for h in fields
    β0 = Dict(SEEDS)[h]
    seconds = @elapsed cert = certify_gap(site_model(h); n0 = N0, β0, rp = site_rp(h),
        blocking = BLOCKING)
    cert["staggered_field"] = string(h)
    cert["seconds"] = seconds
    cert["metadata"] = model_metadata(ROOT, Dict("h" => string(h), "n0" => N0,
        "beta0" => string(β0)); SOURCES..., packages = (CertifiedGaps,), certified = true)
    write_result(certificate_path(h), cert)
    println("CERTIFIED h=$h gap > $(cert["gap_lower_bound_decimal"]) ",
        "($(round(seconds; digits = 1)) s)")
    flush(stdout)
end
# Compact summary of every stored certificate (checked against them by test/runtests.jl).
rows = map(FIELDS) do h
    c = read_result(certificate_path(h), Dict{String,Any})
    Dict("h" => string(h), "n0" => c["seed"]["n0"], "beta0" => c["seed"]["beta0"],
        "eps_s" => c["seed"]["eps_s"], "eps_t" => c["seed"]["eps_t"],
        "gap_lower_bound" => c["gap_lower_bound"],
        "gap_lower_bound_decimal" => c["gap_lower_bound_decimal"],
        "interval_bound_decimal" => c["interval_bound_decimal"],
        "seconds" => c["seconds"])
end
write_result(joinpath(ROOT, "results", "summary.json"), Dict(
    "description" =>
        "Seeds and gap bounds of results/certificate_<h>.json; every even L >= 8 sites",
    "points" => rows))
