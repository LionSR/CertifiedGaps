# Certifies one dimerization delta = d/10: verified spectra of the rings n0 and 2n0 cells
# (complete S^z blocks for small rings, exact-projector symmetry blocks for twelve spins),
# outward log Z, rational seeds and the recurrence. Writes results/certificate_d<d>.json, or
# the path in SPIN1_DIMERIZED_REPLAY_PATH (to replay without overwriting).
#     julia --project=models/spin1_dimerized models/spin1_dimerized/certify.jl 7          (~17 min)
#     julia --project=models/spin1_dimerized models/spin1_dimerized/certify.jl 9 2 7//2   (~1 s)
# Without n0 and beta0 the seed is the best one in results/scan_d<d>.json.
include(joinpath(@__DIR__, "src", "common.jl"))
include(joinpath(@__DIR__, "src", "rp.jl"))
include(joinpath(@__DIR__, "src", "enclosure.jl"))
if abspath(PROGRAM_FILE) == @__FILE__
    certify_selected(parse(Int, ARGS[1]) // 10;
        n0 = length(ARGS) > 1 ? parse(Int, ARGS[2]) : nothing,
        beta0 = length(ARGS) > 2 ? parserat(ARGS[3]) : nothing,
        output_path = get(ENV, "SPIN1_DIMERIZED_REPLAY_PATH", nothing))
end
