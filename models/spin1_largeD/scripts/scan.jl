# Floating-point screening scan (not a certificate): spectra at L = 2..12 for each D, then
# the core's seed search over n0 = 2..6 and beta = 0.125:0.125:20. About 25 s per D.
#
#   julia --project=models/spin1_largeD models/spin1_largeD/scripts/scan.jl [D ...]
#
# Writes results/scan_D*_L*.json (spectra cache) and results/summary_D*.json; both are
# gitignored and summarised by scripts/report.jl. Large blocks use Krylov low spectra with the
# core's cutoff bookkeeping, so nothing here is rigorous.
include(joinpath(@__DIR__, "..", "model.jl"))
using Random
Random.seed!(1701)
root = dirname(@__DIR__)
Ds = isempty(ARGS) ? collect(4:-1//4:1) : [parse(Rational{Int}, x) for x in ARGS]
for D in Ds
    rp_identity(D)
    spectra = Dict{Int,PartialSpectrum}(); records = Any[]
    for L in [2, 3, 4, 5, 6, 8, 10, 12]
        path = joinpath(root, "results", "scan_D$(numerator(D))_$(denominator(D))_L$L.json")
        if isfile(path)
            data = JSON3.read(read(path, String), Dict{String,Any})
        else
            data = compute_spectrum(L, D; low = true)
            data["source_hashes"] = source_hashes(root; dirs = (),
                files = ("model.jl", "scripts/scan.jl"))
            data["metadata"] = run_metadata(Dict("D" => string(D), "L" => L);
                source_sha256 = data["source_hashes"], status = "numerical, not certified")
            write_result(path, data)
        end
        E = Float64.(data["energies"])
        spectra[L] = PartialSpectrum(E, Int(data["dimension"]), data["complete"])
        push!(records, Dict("L" => L, "E0" => E[1], "gap" => E[2]-E[1],
            "seconds" => data["seconds"]))
        println("SCAN D=$D L=$L seconds=$(data["seconds"]) E0=$(E[1]) gap=$(E[2]-E[1])")
        flush(stdout)
    end
    # Odd n are diagnostic only; the spatial moment theorem requires an even seed.
    results = scan_seeds(spectral_defects(spectra), [2, 3, 4, 5, 6], 0.125:0.125:20)
    best = best_seed(filter(x -> iseven(x.n₀), results))
    rows = [Dict("n0" => r.n₀, "beta" => r.β₀, "eps_s" => r.εₛ, "eps_t" => r.εₜ,
        "bound" => r.bound) for r in results]
    diagnostics = [Dict("n0" => n,
        "max_min_purity" => maximum(min(1-r.εₛ, 1-r.εₜ) for r in results if r.n₀ == n),
        "delta" => spectra[2n].E[1] - 2spectra[n].E[1]) for n in [2, 3, 4, 5, 6]]
    data = Dict("D" => string(D), "certified" => false, "rows" => rows,
        "diagnostics" => diagnostics, "sizes" => records,
        "best" => best === nothing ? nothing : Dict("n0" => best.n₀, "beta" => best.β₀,
            "eps_s" => best.εₛ, "eps_t" => best.εₜ, "bound" => best.bound),
        "source_hashes" => source_hashes(root; dirs = (),
            files = ("model.jl", "scripts/scan.jl")))
    write_result(joinpath(root, "results",
        "summary_D$(numerator(D))_$(denominator(D)).json"), data)
    println("BEST D=$D $best"); flush(stdout)
end
