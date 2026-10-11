include(joinpath(@__DIR__, "verify_polynomial_flow.jl"))
using CertifiedGaps.Results: write_result

function reference_refinement(path, output; replay = true)
    replay && verify_flow(path)
    data = JSON3.read(read(path, String))
    initial = JSON3.read(read(String(data.source_file), String))
    seed = PM.scalar_seed(initial)
    kappas, deltas, local_bounds = MN.F[], MN.F[], MN.F[]
    for row in data.runs
        file = joinpath(dirname(path), String(row.witness_file))
        @assert bytes2hex(sha256(read(file))) == row.witness_sha256
        proof = JSON3.read(read(file, String))
        kappa = MN.parse_fraction(proof.normalization)
        push!(kappas, kappa)
        push!(
            deltas,
            MN.parse_fraction(proof.error)+kappa*MN.parse_fraction(
                proof.canonicalization_error_upper,
            ),
        )
        push!(local_bounds, MN.parse_fraction(proof.operator_norm_upper))
    end
    reference = PM.reference_operator_bounds(seed, kappas)
    # Also use the direct Gram difference bound with both independent norms.
    coefficients = local_bounds
    errors = MS.propagate_gram_errors(
        MN.parse_fraction(initial.error_upper),
        kappas,
        deltas,
        vcat(coefficients, reference[end]),
        round_bits = get(data, :error_rounding_bits, nothing),
        reference_bounds = reference,
    )
    result = Dict(
        "source_file"=>path,
        "source_sha256"=>bytes2hex(sha256(read(path))),
        "reference_norm_bounds"=>string.(reference),
        "propagated_error_bounds"=>string.(errors),
        "final_error_display"=>Float64(last(errors)),
        "previous_error_display"=>Float64(
            MN.parse_fraction(last(data.runs).propagated_error_upper),
        ),
        "source_fully_replayed_this_run"=>replay,
        "thermal_trace_certified"=>false,
        "uniform_gap_certified"=>false,
    )
    write_result(output, result)
    println(
        "Analytic reference refinement: $(result["previous_error_display"]) -> $(result["final_error_display"])",
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    reference_refinement(ARGS...)
end
