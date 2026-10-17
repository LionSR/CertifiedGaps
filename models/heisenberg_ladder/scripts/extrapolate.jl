# Heuristic seed-scale estimate at the couplings with no seed through ten rungs (NOT a
# certificate; seconds). Fits δ(n) = δ(4) exp(-(n-4)/ξ) to the even ground-energy mismatches
# δ(n) = E₀(2n) - 2E₀(n) at n = 2, 4 of results/scans/scan_J*.json, models
# Pₛ ≈ exp(-β δ(n)) and P_τ ≈ exp(-6n exp(-β Δ₁₀)) (dilute triplets, Δ₁₀ the ten-rung ED gap),
# and reports the first even n₀ whose modelled seed closes the recurrence on β = 1/4:1/4:32.
# Two sizes cannot validate exponential scaling; no proof uses this.
# Writes results/scans/extrapolation.json.
#
#     julia --project=julia models/heisenberg_ladder/scripts/extrapolate.jl
include(joinpath(@__DIR__, "..", "src", "HeisenbergLadder.jl"))
using .HeisenbergLadder: BASE, HASHED
using CertifiedGaps.Seeds: scan_seeds, best_seed
using CertifiedGaps.Results: write_result, model_metadata
using JSON3

"Estimate the first admissible even seed length for every scan without a seed."
function extrapolate()
    estimates = Any[]
    for tag in ("4p0", "3p5", "3p0", "2p5", "2p0", "1p5")
        s = JSON3.read(read(joinpath(BASE, "results", "scans", "scan_J$(tag).json"), String))
        s.best !== nothing && continue
        ds = Dict(Int(d.n0) => Float64(d.delta) for d in s.diagnostics)
        xi = 2 / log(ds[2] / ds[4]); gap = Float64(s.rings[Symbol("10")].gap)
        prediction = nothing
        for n in 6:2:128
            delta = ds[4] * exp(-(n - 4) / xi)
            defects = (nn, beta) -> (-expm1(-beta * delta), -expm1(-6nn * exp(-beta * gap)))
            b = best_seed(scan_seeds(defects, (n,), 0.25:0.25:32))
            if b !== nothing
                prediction = Dict("n0" => n, "beta" => b.β₀, "eps_s" => b.εₛ, "eps_t" => b.εₜ)
                break
            end
        end
        push!(estimates, Dict("J" => s.J, "xi_from_even_deltas" => xi, "ed_gap_10" => gap,
            "prediction" => prediction))
    end
    write_result(joinpath(BASE, "results", "scans", "extrapolation.json"), Dict(
        "metadata" => model_metadata(BASE, Dict("model" =>
            "delta4 exp(-(n-4)/xi), P_t=exp(-6n exp(-beta gap10))"); HASHED...),
        "estimates" => estimates))
    foreach(println, estimates)
end

extrapolate()
