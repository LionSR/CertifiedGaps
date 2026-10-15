include(joinpath(@__DIR__, "operator_error_probe.jl"))
using CertifiedGaps.Results: write_result

include(joinpath(@__DIR__, "operator_power_screen.jl"))

function refine_operator_error(source, output; powers = 0.0:0.05:0.6)
    @assert abspath(source) != abspath(output)
    text = read(source, String)
    data = JSON3.read(text, Dict{String,Any})
    source_square = data["source_file"]
    square_text = read(source_square, String)
    @assert bytes2hex(sha256(codeunits(square_text))) == data["source_sha256"]
    input, square = MS.decode_witness(JSON3.read(square_text))
    original = OE.decode(data)
    screen = operator_power_screen(input, square, original; powers)
    families = (original.first_gauges, original.second_gauges, original.prefix_gauges)
    changed = ntuple(k -> [s .* g for (s, g) in zip(screen.scales[k], families[k])], 3)
    progress = row -> begin
        if row.site % 16 == 0 || row.site == row.total-1
            println("Refined operator residual: site $(row.site)/$(row.total-1)")
            flush(stdout)
        end
    end
    refined = OE.operator_error_certificate(
        input,
        square;
        first_gauges = changed[1],
        second_gauges = changed[2],
        prefix_gauges = changed[3],
        arithmetic = :interval,
        progress,
    )
    result = OE.encode(refined)
    merge!(
        result,
        Dict(
            "source_file"=>source_square,
            "source_sha256"=>data["source_sha256"],
            "refined_from"=>source,
            "refined_from_sha256"=>bytes2hex(sha256(codeunits(text))),
            "operator_error_upper_display"=>Float64(refined.error),
            "hs_error_upper_display"=>Float64(square.error),
            "selected_gauge_powers"=>screen.powers,
            "numerical_screen_certified"=>false,
            "screen_log_error_before"=>screen.baseline_log,
            "screen_log_error_after"=>screen.best_log,
            "generator_replay_complete"=>false,
            "uniform_gap_certified"=>false,
            "thermal_trace_certified"=>false,
        ),
    )
    write_result(output, result)
    @assert OE.verify_operator_error(
        input,
        square,
        refined;
        residual_arithmetic = :interval,
        progress,
    )
    result["generator_replay_complete"] = true
    write_result(output, result)
    println(
        "Refined operator error <= $(Float64(refined.error)); previous <= $(Float64(original.error))",
    )
    return result
end

if abspath(PROGRAM_FILE) == @__FILE__
    @assert length(ARGS) == 2
    refine_operator_error(ARGS...)
end
