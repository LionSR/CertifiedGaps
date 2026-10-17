# Floating-point screening scan (not a certificate): spectra at L = 2..12 for each D from the
# core's charge-sector exact diagonalization (dense up to DENSE_MAX states per sector, Lanczos
# lowest KEEP levels above that, with the core's cutoff bookkeeping), then the core's seed search
# over n0 = 2..6 and beta = 0.125:0.125:20.
#
#   julia --project=julia models/spin1_largeD/scripts/scan.jl [D ...]
#
# Writes results/scan_D*_L*.json (spectra cache) and results/summary_D*.json; both are
# gitignored and summarized by scripts/report.jl. Nothing here is rigorous.
using CertifiedGaps, CertifiedGaps.ExactDiagonalization, CertifiedGaps.Seeds
using CertifiedGaps.Results, JSON3, LinearAlgebra, Random
include(joinpath(dirname(@__DIR__), "src", "model.jl"))
using .LargeDModel
BLAS.set_num_threads(1)
Random.seed!(1701)
const ROOT = dirname(@__DIR__)
const DENSE_MAX, KEEP = 4000, 400

"The ring of `ring_model(D)` as a core floating-point `Chain` (same bond, on-site term, charge)."
function chain(D)
    m = ring_model(D)
    return Chain(3, ring_bonds(Float64.(real(m.bond))); onsite = Float64.(real(m.onsite)),
        charge = m.charge)
end

"Spectrum of the ring of L sites at D, cached in results/scan_D*_L*.json."
function cached_spectrum(D, L)
    path = joinpath(ROOT, "results", "scan_$(stem(D))_L$L.json")
    if isfile(path)
        data = JSON3.read(read(path, String), Dict{String,Any})
        return PartialSpectrum(Float64.(data["energies"]), Int(data["dimension"]),
            data["complete"]), data["seconds"]
    end
    seconds = @elapsed s = spectrum(chain(D), L; spin_flip = true, dense_max = DENSE_MAX,
        k = KEEP)
    (; source_sha256) = model_provenance(ROOT; dirs = (),
        files = ("src/model.jl", "scripts/scan.jl"))
    write_result(path, Dict("D" => string(D), "L" => L, "energies" => s.E,
        "dimension" => s.N, "complete" => s.exact, "seconds" => seconds,
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
    rows = [Dict("n0" => r.n₀, "beta" => r.β₀, "eps_s" => r.εₛ, "eps_t" => r.εₜ,
        "bound" => r.bound) for r in results]
    diagnostics = [Dict("n0" => n,
        "max_min_purity" => maximum(min(1 - r.εₛ, 1 - r.εₜ) for r in results if r.n₀ == n),
        "delta" => spectra[2n].E[1] - 2spectra[n].E[1]) for n in [2, 3, 4, 5, 6]]
    data = Dict("D" => string(D), "certified" => false, "rows" => rows,
        "diagnostics" => diagnostics, "sizes" => records,
        "best" => best === nothing ? nothing : Dict("n0" => best.n₀, "beta" => best.β₀,
            "eps_s" => best.εₛ, "eps_t" => best.εₜ, "bound" => best.bound))
    write_result(joinpath(ROOT, "results", "summary_$(stem(D)).json"), data)
    println("BEST D=$D $best"); flush(stdout)
end
