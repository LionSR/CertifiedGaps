# A different rigorous upper-bound recurrence. If ordered eigenvalues of M
# dominate those of A >= 0, and C + delta*I >= M^2, then the new ordered
# eigenvalues dominate those of A^2. No matrix-monotonicity of squaring is used.
include(joinpath(@__DIR__, "repeated_l2_probe.jl"))
using CertifiedGaps.Results: write_result

function identity_amplitude()
    result = zeros(F, 625, 625)
    for r in 0:624
        result[location(r, r)...] = 1
    end
    return result
end

function run_loewner(output; rankcaps = (16, 32, 64), beta_target = 4, doublings = 5)
    result = Dict{String,Any}(
        "status"=>"four-site-certified-spectral-majorant-compression-probe",
        "complete"=>false,
        "uniform_gap_certified"=>false,
        "length"=>4,
        "beta_target"=>beta_target,
        "doublings"=>doublings,
        "steps_final"=>2^(doublings+1),
        "runs"=>Any[],
        "runner_sha256"=>bytes2hex(sha256(read(@__FILE__))),
        "source_sha256"=>Dict(
            name=>bytes2hex(sha256(read(joinpath(@__DIR__, name)))) for
            name in ("repeated_l2_probe.jl", "positive_compression_probe.jl")
        ),
    )
    references = [amplitude_interval(beta_target//2^(doublings-k), 2^(k+1)) for k in 0:doublings]
    identity = identity_amplitude()
    for rankcap in rankcaps
        lower, upper, _ = first(references)
        rounding = sqrt_upper(norm2((upper-lower)/2))
        M = hermitize((lower+upper)/2) + rounding*identity
        started = time()
        for k in 0:doublings
            if k > 0
                lower, upper = squared_interval(M)
                M = (lower+upper)/2
                rounding = sqrt_upper(norm2((upper-lower)/2))
            else
                rounding = zero(F)
            end
            C, rank = svd_candidate(M, rankcap)
            delta = operator_norm_upper(M-C)+rounding
            grid = big(2)^80
            delta = cld(numerator(delta)*grid, denominator(delta))//grid
            M = C+delta*identity
            trace_upper = norm2(M)
            reference = last(references[k+1])
            @assert reference < trace_upper
            row = Dict(
                "rankcap_per_charge"=>rankcap,
                "beta"=>string(beta_target//2^(doublings-k)),
                "doublings_completed"=>k,
                "represented_rank_upper"=>min(625, rank+1),
                "identity_correction"=>string(delta),
                "identity_correction_display"=>Float64(delta),
                "thermal_trace_upper"=>string(trace_upper),
                "direct_trotter_trace_ball"=>string(reference),
                "cumulative_inflation_display"=>Float64(trace_upper)/Float64(reference),
                "seconds"=>time()-started,
            )
            push!(result["runs"], row)
            write_result(output, result)
            println(
                "spectral compression cap=$rankcap stage=$k beta=$(row["beta"]) upper/reference=$(row["cumulative_inflation_display"]) identity=$(Float64(delta)) seconds=$(time()-started)",
            )
            flush(stdout)
        end
    end
    result["complete"] = true
    result["completed_utc"] = string(now(UTC))
    write_result(output, result)
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_loewner(
        isempty(ARGS) ? joinpath(@__DIR__, "../../results/loewner-compression.json") : only(ARGS),
    )
end
