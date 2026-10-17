# Numerical seed scan (floating point, not certified): spectra of physical-frame rings up to maxL
# sites from the core exact diagonalization (Sᶻ sectors; Lanczos lowest 80 levels above
# dimension 2500), cached in results/spectrum_<h>_<L>.json, and a β grid written to
# results/scan_<h>_<maxL>.json (both gitignored). The certified points are in src/common.jl.
#     julia --project=julia models/heisenberg_staggered_field/scripts/scan.jl 16 [h]
include(joinpath(dirname(@__DIR__), "src", "common.jl"))
using CertifiedGaps.ExactDiagonalization, CertifiedGaps.Seeds, JSON3
maxL = length(ARGS) > 0 ? parse(Int, ARGS[1]) : 16
fields = length(ARGS) > 1 ? [parse_rational(ARGS[2])] : FIELDS
for h in fields
    chain = CertifiedGaps.Models.staggered_heisenberg(h).chain
    spectra = Dict{Int,PartialSpectrum}()
    timings = Dict{Int,Float64}()
    lengths = sort!(unique(vcat(collect(2:2:maxL÷2), collect(4:4:maxL))))
    for L in lengths
        path = joinpath(ROOT, "results", "spectrum_$(tag(h))_$(L).json")
        if isfile(path)
            d = JSON3.read(read(path, String), Dict{String,Any})
            spectra[L] = PartialSpectrum(Float64.(d["eigenvalues"]), 2^L, d["complete"])
            timings[L] = d["seconds"]
        else
            timings[L] = @elapsed spectra[L] = spectrum(chain, L; dense_max = 2500, k = 80)
            write_result(path, Dict("h" => string(h), "L" => L, "eigenvalues" => spectra[L].E,
                "complete" => spectra[L].exact, "seconds" => timings[L]))
        end
    end
    rows = scan_seeds(spectral_defects(spectra), 2:2:maxL÷2, 0.25:0.25:32)
    best = best_seed(rows)
    table = [Dict("n0" => n,
        "max_min_purity" => maximum(min(1 - r.εₛ, 1 - r.εₜ) for r in rows if r.n₀ == n),
        "delta" => spectra[2n].E[1] - 2spectra[n].E[1]) for n in 2:2:maxL÷2]
    bestdict = best === nothing ? nothing : Dict("n0" => best.n₀, "beta" => best.β₀,
        "epsilon_s" => best.εₛ, "epsilon_t" => best.εₜ, "bound" => best.bound)
    top = spectra[maximum(lengths)].E
    write_result(joinpath(ROOT, "results", "scan_$(tag(h))_$(maxL).json"), Dict(
        "status" => "floating-point scan, not certified", "h" => string(h),
        "best" => bestdict, "per_n" => table, "timings" => timings,
        "reference_gap" => top[2] - top[1],
        "grid" => [Dict("n0" => r.n₀, "beta" => r.β₀, "es" => r.εₛ, "et" => r.εₜ,
            "bound" => r.bound) for r in rows]))
    println("BEST h=$h $bestdict")
    flush(stdout)
end
