# β₀ scan of the coupled seeds at n₀ = 6 cells (floating-point screening of rigorous inputs):
# for each δ of RITZ_SEEDS, the coupled defects from the complete verified 12-spin spectrum and
# the stored 24-spin Ritz bounds (results/ritz_bounds_12_cells_delta_*.json.gz) on the grid
# β₀ ∈ [lo, hi] in steps of 1/8, and the interval bound `certified_gap_bound`. RITZ_SEEDS holds
# the maximizers. Also reports the same scan with the bounds truncated to the 8 lowest levels per
# block, as a check that the omitted weight is negligible.
#
#   julia --project=code/julia code/models/dimerized_heisenberg/scripts/ritz_scan.jl [p q lo hi]
using CertifiedGaps, CertifiedGaps.VerifiedSpectra, CertifiedGaps.RitzBounds, CertifiedGaps.Seeds
using Printf
include(joinpath(@__DIR__, "..", "src", "DimerizedHeisenberg.jl"))
include(joinpath(@__DIR__, "..", "src", "RitzSeeds.jl"))
using .DimerizedHeisenberg, .RitzSeeds

slug(x) = denominator(x) == 1 ? "$(numerator(x))" : "$(numerator(x))_$(denominator(x))"

function scan(δ, lo, hi)
    large = read_partial_spectrum(joinpath(@__DIR__, "..", "results",
        "ritz_bounds_12_cells_delta_$(slug(δ)).json.gz"))
    short = [PartialBlock(b.upper[1:min(8, end)], b.multiplicity, b.label) for b in large]
    small = verified_ring_spectrum(rotated_model(δ), RITZ_N0)
    best = (0.0, 0 // 1)
    for k in lo:hi
        β = big(k) // 8
        zn = log_partition_enclosure(small, β)
        row = map((large, short)) do bounds
            c = coupled_seed_defects(zn, log_partition_lower(bounds, β),
                log_partition_lower(bounds, 2β))
            (c, certified_gap_bound(c.eps_s, c.eps_t, β).bound)
        end
        (c, b), (_, b8) = row
        @printf("delta=%s beta0=%s eps_s=%.6f eps_t=%.6f bound=%s (8 levels: %s)\n", δ, β,
            Float64(c.eps_s), Float64(c.eps_t), something(b, "-"), something(b8, "-"))
        b !== nothing && b > best[1] && (best = (b, β))
    end
    println("best delta=$δ: beta0=$(best[2]) bound=$(best[1])")
end

if isempty(ARGS)
    for (δ, β0) in sort!(collect(RITZ_SEEDS); rev = true)
        k = Int(8β0)
        scan(δ, k - 8, k + 8)
    end
else
    scan(parse(BigInt, ARGS[1]) // parse(BigInt, ARGS[2]), parse(Int, ARGS[3]),
        parse(Int, ARGS[4]))
end
