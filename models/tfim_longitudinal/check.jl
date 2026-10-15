# Checker for the 30 longitudinal-field TFIM certificates in results/certificate_g*_h*.json.gz.
#
#     julia --project=models/tfim_longitudinal models/tfim_longitudinal/check.jl [--full]
#
# Fast mode (default, ~20 s): exact RP identity, momentum-sector dimensions, complete spectra,
# the interval partition sums and rational seeds recomputed from the stored eigenvalue
# intervals, and the recurrence bound. It trusts the stored eigenvalue radii.
# `--full` (~6 min): also rebuilds every block and its rigorous enclosure with the current code
# and checks, in exact rationals, that each fresh interval lies inside the stored one, so the
# radii are re-proved rather than trusted. Any failure throws. The checks themselves are
# `check_result` and `check_grid` in `src/checks.jl`, which the tests reuse.
include(joinpath(@__DIR__, "src", "LongitudinalCertificate.jl"))
using .LongitudinalCertificate
const C=LongitudinalCertificate

full="--full" in ARGS
paths=C.result_paths()
C.check_grid(paths)
for path in paths
    println(basename(path),": ",C.check_result(path; full))
end
println("All $(length(paths)) certificates verified", full ? " (full replay)." : ".")
