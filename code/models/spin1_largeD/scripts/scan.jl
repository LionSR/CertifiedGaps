# Floating-point screening scan (not a certificate): complete spectra at L = 2..12 for each D
# from the core's unverified block eigenvalues (VerifiedSpectra.float_ring_spectrum, the same
# charge, momentum, flip and reflection blocks as the certificates), then the core's seed search
# over n0 = 2..6 and beta = 0.125:0.125:20.
#
#   julia --project=julia models/spin1_largeD/scripts/scan.jl [D ...]
#
# Writes results/scan_D*_L*.json (spectra cache) and results/summary_D*.json; both are
# gitignored and summarized by scripts/report.jl. Nothing here is rigorous.
using CertifiedGaps.ExactDiagonalization: PartialSpectrum
using CertifiedGaps.VerifiedSpectra: float_ring_spectrum
using CertifiedGaps.Seeds: scan_seeds, spectral_defects, best_seed, seed_record, purity_peaks
using CertifiedGaps.Results, JSON3, LinearAlgebra
include(joinpath(dirname(@__DIR__), "src", "model.jl"))
using .LargeDModel
BLAS.set_num_threads(1)
const ROOT = dirname(@__DIR__)

"Spectrum of the ring of L sites at D, cached in results/scan_D*_L*.json."
function cached_spectrum(D, L)
    path = joinpath(ROOT, "results", "scan_$(stem(D))_L$L.json")
    if isfile(path)
        data = JSON3.read(read(path, String), Dict{String,Any})
        return PartialSpectrum(Float64.(data["energies"]), Int(data["dimension"])),
        data["seconds"]
    end
    seconds = @elapsed s = float_ring_spectrum(ring_model(D), L)
    (; source_sha256) = model_provenance(ROOT; dirs = (),
        files = ("src/model.jl", "scripts/scan.jl"))
    write_result(path, Dict("D" => string(D), "L" => L, "energies" => s.E,
        "dimension" => s.N, "seconds" => seconds,
        "metadata" => run_metadata(Dict("D" => string(D), "L" => L); source_sha256,
            status = "numerical, not certified")))
    return s, seconds
end

Ds = isempty(ARGS) ? collect(4:-1//4:1) : [parse(Rational{Int}, x) for x in ARGS]
for D in Ds
    spectra = Dict{Int,PartialSpectrum}(); records = Any[]
    for L in [2, 3, 4, 5, 6, 8, 10, 12]
        s, seconds = cached_spectrum(D, L)
        spectra[L] = s
        push!(records, Dict("L" => L, "E0" => s.E[1], "gap" => s.E[2] - s.E[1],
            "seconds" => seconds))
        println("SCAN D=$D L=$L seconds=$seconds E0=$(s.E[1]) gap=$(s.E[2]-s.E[1])")
        flush(stdout)
    end
    # Odd n are diagnostic only; the spatial moment theorem requires an even seed.
    results = scan_seeds(spectral_defects(spectra), [2, 3, 4, 5, 6], 0.125:0.125:20)
    best = best_seed(filter(x -> iseven(x.n₀), results))
    diagnostics = map(purity_peaks(results)) do peak
        n = peak["n0"]
        merge(peak, Dict("delta" => spectra[2n].E[1] - 2spectra[n].E[1]))
    end
    data = Dict("D" => string(D), "certified" => false, "rows" => seed_record.(results),
        "diagnostics" => diagnostics, "sizes" => records,
        "best" => best === nothing ? nothing : seed_record(best))
    write_result(joinpath(ROOT, "results", "summary_$(stem(D)).json"), data)
    println("BEST D=$D $best"); flush(stdout)
end
