# Stage 2: write results/certificate.json and results/finite_bridge.json (about 2 min).
#
#   OPENBLAS_NUM_THREADS=1 julia --project=julia \
#       models/spin1_largeD/scripts/transfer_certify.jl
#
# Reads only the stored, hash-checked L = 6, 12 verified spectra written by certify.jl and
# recomputes the small L = 8, 10 bridge spectra. No floating-point shortcut decides acceptance:
# every candidate is evaluated with interval log Z and exact rational defects, and a seed is kept
# only if `certified_gap_bound` accepts it.
include(joinpath(dirname(@__DIR__), "src", "TransferCertificate.jl"))
using .TransferCertificate, CertifiedGaps, CertifiedGaps.Results, IntervalArithmetic
const TC = TransferCertificate

"Best successful row over cached log Z enclosures for the interval [lo, hi], or `nothing`."
function best_transported(cache, Dref, lo, hi; max_lifts = 5)
    best = nothing
    for (beta, base) in cache
        for row in TC.evaluate_logs(base, beta, Dref, lo, hi; max_lifts)
            row.result.bound === nothing && continue
            if best === nothing || row.result.bound > best.result.bound
                best = row
            end
        end
    end
    return best
end

"""
    certify_transfer()

Re-certify the stage-1 seeds, search spatial lifts at each reference D, certify the fixed
transported D intervals `TC.TRANSPORT_LO`, and record by bisection how far below each reference
D the transport would still succeed (diagnostic); write `results/certificate.json`.
"""
function certify_transfer()
    started = time()
    records, diagnostics = Any[], Any[]
    for (D, beta0) in TC.POINTS
        a, b = TC.reference_spectrum(D, 6), TC.reference_spectrum(D, 12)
        direct = only(TC.evaluate(a, b, beta0; max_lifts = 0))
        push!(records, TC.record(direct, b, "reference-D$D"))
        if D == 4
            # The survey's quoted seed (n0, beta0) = (6, 3), now with an interval bound.
            survey = only(TC.evaluate(a, b, 3 // 1; max_lifts = 0))
            push!(records, TC.record(survey, b, "survey-D4-beta3"))
        end
        cache = [(beta, [TC.logz(a, beta), TC.logz(b, beta),
            TC.logz(a, 2beta), TC.logz(b, 2beta)]) for beta in
                 sort(unique(vcat(collect(3//2:1//40:6), [beta0])))]
        for lifts in 0:5
            best = nothing
            for (beta, base) in cache
                row = TC.evaluate_logs(base, beta, D, D, D; max_lifts = lifts)[end]
                row.result.bound === nothing && continue
                if best === nothing || row.result.bound > best.result.bound
                    best = row
                end
            end
            push!(diagnostics, Dict("D" => string(D), "lifts" => string(lifts),
                "n0" => string(6 * 2^lifts), "best" => best === nothing ? nothing :
                    TC.record(best, b, "space-lift-$lifts-D$D")))
            best === nothing || lifts == 0 || push!(records,
                TC.record(best, b, "space-lift-$lifts-D$D"))
        end
        best = best_transported(cache, D, TC.TRANSPORT_LO[D], D)
        best === nothing && error("the fixed transported interval no longer certifies")
        push!(records, TC.record(best, b, "transport-D$D"))
        # The Weyl radius grows monotonically with |D - Dref|, so integer bisection on a
        # 1e-7 grid finds the last grid point at which some cached seed still succeeds.
        grid = 10^7
        success, failure = 0, grid ÷ 10
        best_transported(cache, D, D - failure // grid, D) === nothing ||
            error("increase the search bracket")
        while failure - success > 1
            middle = (success + failure) ÷ 2
            found = best_transported(cache, D, D - middle // grid, D)
            found === nothing ? (failure = middle) : (success = middle)
        end
        push!(diagnostics, Dict("D_reference" => string(D),
            "certified_interval_lower" => string(TC.TRANSPORT_LO[D]),
            "smallest_successful_transport_grid_D" => string(D - success // grid),
            "next_transport_grid_D" => string(D - failure // grid),
            "D_grid_step" => string(1 // grid),
            "beta_grid" => "3/2:1/40:6, plus the stage-1 beta0",
            "interpretation" => "reach of these Weyl/transfer upper bounds only"))
        println("Dref=$D transported to [$(best.lo), $D]: n=$(best.n0), beta=$(best.beta), " *
                "bound=$(best.result.bound); bisection reaches $(D - success // grid)")
        flush(stdout)
    end
    data = Dict("format" => TC.FORMAT, "cases" => records,
        "search_diagnostics" => diagnostics, "input_hashes" => TC.input_hashes(),
        "wall_seconds" => time() - started,
        "metadata" => TC.metadata(Dict("largest_ED_ring" => "12")))
    write_result(joinpath(TC.ROOT, "results", "certificate.json"), data)
    return data
end

if abspath(PROGRAM_FILE) == @__FILE__
    certify_transfer()
    bridge = TC.bridge_certificate()
    write_result(joinpath(TC.ROOT, "results", "finite_bridge.json"), bridge)
    println("finite bridge: ", [(r["L"], r["gap_lower"]) for r in bridge["rings"]])
end
