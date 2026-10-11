# Reference experiment for repeated certified compression, not a large-chain
# algorithm. The dense four-site operator is resolved into charge sectors.
include(joinpath(@__DIR__, "positive_compression_probe.jl"))
using CertifiedGaps.Results: write_result

function square_majorant(M; bits = 80)
    RR = Nemo.ArbField(192)
    result = zeros(F, 625, 625)
    powers = [5^(j-1) for j in 1:4]
    for mag in -8:8
        states = [code for code in 0:624 if sum(mod(div(code, p), 5)-2 for p in powers) == mag]
        block =
            [M[mod(r, 25)+25mod(c, 25)+1, div(r, 25)+25div(c, 25)+1] for r in states, c in states]
        squared = Nemo.matrix(RR, block)^2
        for (i, r) in enumerate(states), (j, c) in enumerate(states)
            _, upper = PB.rational_bracket(squared[i, j]; bits)
            result[mod(r, 25)+25mod(c, 25)+1, div(r, 25)+25div(c, 25)+1] = max(upper, zero(F))
        end
    end
    return result
end

function compress_majorant(M, rankcap; iterations = 3000, underweight = 4096.0)
    metric = [local_metric(i) for i in 0:624]
    groups = Dict(q=>findall(i->local_charge(i-1)==q, 1:625) for q in -8:8)
    bound = zeros(F, 625, 625)
    totalrank = 0
    for q in -8:8
        rows, cols = groups[q], groups[-q]
        block = M[rows, cols]
        iszero(block) && continue
        c = PB.upper_factorization(
            block,
            min(rankcap, length(rows), length(cols));
            iterations,
            underweight,
            row_metric = metric[rows],
            col_metric = metric[cols],
        )
        bound[rows, cols] = c.bound
        totalrank += c.rank
    end
    @assert all(bound .>= M)
    return bound, totalrank
end

function weighted_norm2(M)
    metric = [local_metric(i) for i in 0:624]
    return sum(M[i, j]^2*metric[i]*metric[j] for i in 1:625, j in 1:625)
end

function run_repeated(
    output = joinpath(@__DIR__, "../results/repeated-positive-compression.json");
    rankcaps = (32, 64),
    beta_target = 4,
    doublings = 5,
)
    result = Dict{String,Any}(
        "status"=>"four-site-repeated-positive-compression-probe",
        "complete"=>false,
        "uniform_gap_certified"=>false,
        "length"=>4,
        "beta_target"=>beta_target,
        "doublings"=>doublings,
        "steps_final"=>2^(doublings+1),
        "runs"=>Any[],
        "runner_sha256"=>bytes2hex(sha256(read(@__FILE__))),
        "source_sha256"=>Dict(
            name=>bytes2hex(sha256(read(joinpath(@__DIR__, "../src", name)))) for
            name in ("../../../julia/src/CertifiedMPO/PositiveBounds.jl", "ThermalBounds.jl")
        ),
        "reference_runner_sha256"=>bytes2hex(
            sha256(read(joinpath(@__DIR__, "positive_compression_probe.jl"))),
        ),
    )
    references = [amplitude_upper(beta_target//2^(doublings-k), 2^(k+1)) for k in 0:doublings]
    for rankcap in rankcaps
        M = first(references[1])
        started = time()
        for k in 0:doublings
            k > 0 && (M = square_majorant(M))
            before = weighted_norm2(M)
            M, rank = compress_majorant(M, rankcap)
            upper = weighted_norm2(M)
            # Both are independently rounded from the same exact amplitude.
            # The norm check compares against the direct uncompressed power.
            reference = last(references[k+1])
            @assert reference < upper
            row = Dict(
                "rankcap_per_charge"=>rankcap,
                "doublings_completed"=>k,
                "beta"=>string(beta_target//2^(doublings-k)),
                "bond_dimension"=>rank,
                "certified_entrywise_majorant"=>true,
                "thermal_trace_upper_rational"=>string(upper),
                "direct_trotter_trace_ball"=>string(reference),
                "cumulative_inflation_display"=>Float64(upper)/Float64(reference),
                "this_compression_inflation_display"=>Float64(upper/before),
                "seconds"=>time()-started,
            )
            push!(result["runs"], row)
            write_result(output, result)
            println(
                "repeated positive cap=$rankcap stage=$k beta=$(row["beta"]) D=$rank upper/reference=$(row["cumulative_inflation_display"]) elapsed=$(time()-started)s",
            )
            flush(stdout)
        end
    end
    result["complete"] = true
    result["completed_utc"] = string(now(UTC))
    write_result(output, result)
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_repeated(
        isempty(ARGS) ? joinpath(@__DIR__, "../results/repeated-positive-compression.json") :
        only(ARGS),
    )
end
