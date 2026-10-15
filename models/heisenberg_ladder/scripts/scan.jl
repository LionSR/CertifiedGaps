# Floating-point seed scan at one coupling (NOT a certificate). Dense complete spectra of the
# rings of 2–6, 8 and 10 rungs, then every seed n₀ = 2:5, β₀ = 1/4:1/4:16. The ten-rung ring
# (2^20 states, largest block 8398) takes 20–35 minutes; peak memory is a few GiB per process.
#
#     OPENBLAS_NUM_THREADS=1 julia --project=models/heisenberg_ladder \
#         models/heisenberg_ladder/scripts/scan.jl 7//2
#
# Writes results/scans/spectrum_J<tag>_L<L>.json (sector table, log Z on β = 1/4:1/4:32 and
# the lowest 32 levels; the full level list is not stored) and results/scans/scan_J<tag>.json.
# Odd n₀ are frustrated diagnostics only: the criterion needs even n₀.
include(joinpath(@__DIR__, "..", "src", "HeisenbergLadder.jl"))
using .HeisenbergLadder: BASE, full_spectrum, provenance
using CertifiedGaps.ExactDiagonalization: PartialSpectrum, log_partition
using CertifiedGaps.Seeds: scan_seeds, spectral_defects, best_seed
using CertifiedGaps.Results: write_result

"File tag of a coupling, e.g. 7//2 -> \"3p5\"."
tag(J) = replace(string(Float64(J)), "." => "p")

seedrow(r) = r === nothing ? nothing : Dict("n0" => r.n₀, "beta" => r.β₀, "eps_s" => r.εₛ,
    "eps_t" => r.εₜ, "bound" => r.bound, "eligible_even_seed" => iseven(r.n₀))

"Run the scan at coupling `J` and write its spectrum and scan files."
function scan(J; rings = (2, 3, 4, 5, 6, 8, 10), ns = 2:5, out = joinpath(BASE, "results", "scans"))
    start = time(); spectra = Dict{Int,PartialSpectrum}(); records = Dict{String,Any}()
    for L in rings
        s = full_spectrum(L, J); es = s.energies
        write_result(joinpath(out, "spectrum_J$(tag(J))_L$(L).json"), Dict(
            "metadata" => provenance(Dict("J" => string(J), "L" => L)),
            "J" => string(J), "L" => L, "dimension" => s.dimension, "sectors" => s.rows,
            "seconds" => s.seconds, "gap" => es[2] - es[1],
            "lowest_32_numerical_levels" => es[1:min(32, length(es))],
            "logZ_grid" => [Dict("beta" => string(Rational{Int}(beta)),
                "logZ" => log_partition(es, beta)) for beta in 0.25:0.25:32]))
        spectra[L] = PartialSpectrum(es, 4^L)
        records[string(L)] = Dict("ground_energy" => es[1], "gap" => es[2] - es[1],
            "dimension" => s.dimension, "seconds" => s.seconds,
            "maximum_block_dimension" => maximum(r["dimension"] for r in s.rows))
    end
    rows = scan_seeds(spectral_defects(spectra), ns, 0.25:0.25:16.0)
    eligible = filter(r -> iseven(r.n₀), rows)
    diagnostics = Any[]
    for n in ns
        rs = filter(r -> r.n₀ == n, rows); best = argmax(r -> min(1 - r.εₛ, 1 - r.εₜ), rs)
        push!(diagnostics, Dict("n0" => n, "max_min_purity" => min(1 - best.εₛ, 1 - best.εₜ),
            "beta_at_max_min" => best.β₀, "delta" => spectra[2n].E[1] - 2spectra[n].E[1],
            "best_bound_seed" => seedrow(best_seed(rs)), "eligible_even_seed" => iseven(n)))
    end
    write_result(joinpath(out, "scan_J$(tag(J)).json"), Dict(
        "metadata" => provenance(Dict("J" => string(J), "beta_grid" => "1//4:1//4:16//1")),
        "J" => string(J), "best" => seedrow(best_seed(eligible)), "rings" => records,
        "diagnostics" => diagnostics, "grid" => seedrow.(rows), "seconds" => time() - start))
    println("BEST ", best_seed(eligible))
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && scan(parse(Rational{Int}, ARGS[1]))
