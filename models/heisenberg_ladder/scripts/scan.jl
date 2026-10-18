# Floating-point seed scan at one coupling (NOT a certificate). Unverified dense eigenvalues of
# the core symmetry blocks of `ladder_model` (`VerifiedSpectra.float_ring_spectrum`) of the
# rings of 2–6, 8 and 10 rungs, then every seed n₀ = 2:5, β₀ = 1/4:1/4:16.
#
#     OPENBLAS_NUM_THREADS=1 julia --project=julia \
#         models/heisenberg_ladder/scripts/scan.jl 7//2
#
# Writes results/scans/spectrum_J<tag>_L<L>.json (log Z on β = 1/4:1/4:32 and the lowest 32
# levels; the full level list is not stored) and results/scans/scan_J<tag>.json.
# Odd n₀ are frustrated diagnostics only: the criterion needs even n₀.
include(joinpath(@__DIR__, "..", "src", "HeisenbergLadder.jl"))
using .HeisenbergLadder: BASE, ladder_model
using CertifiedGaps.VerifiedSpectra: float_ring_spectrum
using CertifiedGaps.ExactDiagonalization: PartialSpectrum, log_partition
using CertifiedGaps.Seeds: scan_seeds, spectral_defects, best_seed, seed_record, purity_peaks
using CertifiedGaps.Results: write_result, model_metadata

"File tag of a coupling, e.g. 7//2 -> \"3p5\"."
tag(J) = replace(string(Float64(J)), "." => "p")

"Run the scan at coupling `J` and write its spectrum and scan files."
function scan(J; rings = (2, 3, 4, 5, 6, 8, 10), ns = 2:5,
    out = joinpath(BASE, "results", "scans"))
    start = time(); spectra = Dict{Int,PartialSpectrum}(); records = Dict{String,Any}()
    model = ladder_model(J)
    for L in rings
        seconds = @elapsed s = float_ring_spectrum(model, L)
        s.exact || error("the blocks do not cover the ring of $L rungs")
        es = s.E
        write_result(joinpath(out, "spectrum_J$(tag(J))_L$(L).json"), Dict(
            "metadata" => model_metadata(BASE, Dict("J" => string(J), "L" => L)),
            "J" => string(J), "L" => L, "dimension" => s.N, "seconds" => seconds,
            "gap" => es[2] - es[1], "lowest_32_numerical_levels" => es[1:min(32, length(es))],
            "logZ_grid" => [Dict("beta" => string(Rational{Int}(beta)),
                "logZ" => log_partition(es, beta)) for beta in 0.25:0.25:32]))
        spectra[L] = s
        records[string(L)] = Dict("ground_energy" => es[1], "gap" => es[2] - es[1],
            "dimension" => s.N, "seconds" => seconds)
    end
    rows = scan_seeds(spectral_defects(spectra), ns, 0.25:0.25:16.0)
    best = best_seed(filter(r -> iseven(r.n₀), rows))
    # The ground-energy mismatch δ(n) = E₀(2n) − 2E₀(n) feeds scripts/extrapolate.jl.
    diagnostics = [merge(p, Dict("delta" => spectra[2p["n0"]].E[1] - 2spectra[p["n0"]].E[1],
        "eligible_even_seed" => iseven(p["n0"]))) for p in purity_peaks(rows)]
    write_result(joinpath(out, "scan_J$(tag(J)).json"), Dict(
        "metadata" => model_metadata(BASE,
            Dict("J" => string(J), "beta_grid" => "1//4:1//4:16//1")),
        "J" => string(J), "best" => best === nothing ? nothing : seed_record(best),
        "rings" => records, "diagnostics" => diagnostics, "grid" => seed_record.(rows),
        "seconds" => time() - start))
    println("BEST ", best)
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && scan(parse(Rational{Int}, ARGS[1]))
