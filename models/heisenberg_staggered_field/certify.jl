# Certifies the staggered-field points with the shared core: for each h in `SEEDS`
# (src/common.jl), `GapCertificates.certify_point` on the one-site rotated-frame model at
# (n₀, β₀) = (8, β₀), with the seed spectra split by the two-site U(1) cell charge `BLOCKING`.
# Writes results/certificate_<h>.json and results/summary.json.
#     julia --project=julia models/heisenberg_staggered_field/certify.jl [h]   (default: all five)
include(joinpath(@__DIR__, "src", "common.jl"))
fields = isempty(ARGS) ? FIELDS : [parse_rational(ARGS[1])]
for h in fields
    seconds = @elapsed cert = certify_point(certificate_path(h), site_model(h); point(h)...)
    println("CERTIFIED h=$h gap > $(cert["gap_lower_bound_decimal"]) ",
        "($(round(seconds; digits = 1)) s)")
    flush(stdout)
end
# Compact summary of every stored certificate (check.jl compares it with them).
write_result(joinpath(ROOT, "results", "summary.json"), Dict(
    "description" =>
        "Seeds and gap bounds of results/certificate_<h>.json; every even L >= 8 sites",
    "points" => [certificate_row(read_result(certificate_path(h), Dict{String,Any}))
                 for h in FIELDS]))
