include(joinpath(@__DIR__, "polynomial_seed_probe.jl"))
import ThermalMPOGap
using ThermalMPOGap.Results: write_result
using ThermalMPOGap.CertifiedMPO.GaugeSearch
const GS = ThermalMPOGap.CertifiedMPO.GaugeSearch
import Nemo

function flow_probe(
    source,
    output;
    steps = 2,
    output_dimension = 32,
    retained_dimension = 16,
    gauge_iterations = 100,
    gram = true,
    reflect_output = true,
    shift = MN.F(47612463, 10000000),
    error_round_bits = 128,
    gauge_mode = :diagonal,
    tail_blocksize = 2,
    candidate_arithmetic = :interval,
    error_method = :triangle,
    compression_bits = 48,
    project_bridge = false,
    conserve_charge = false,
    resume_from = nothing,
    reuse_square_from = nothing,
    reuse_gauge_power = nothing,
    gauge_refiner = nothing,
)
    physical_charges = conserve_charge ? collect(-2:2) : nothing
    conserve_charge &&
        !isnothing(reuse_square_from) &&
        throw(ArgumentError("Charge-preserving generation cannot reuse a fixed dense compression"))
    conserve_charge &&
        !reflect_output &&
        throw(ArgumentError("Charge-preserving flow currently requires reflect_output=true"))
    data = JSON3.read(read(source, String))
    n, power = Int(data.length), Int(data.power)
    seed = PM.scalar_seed(data)
    @assert seed.remainder == MN.parse_fraction(data.remainder)
    B = MN.decode_cores(data.candidate)
    epsilon = MN.parse_fraction(data.error_upper)
    initial_epsilon = epsilon
    println("Checking the saved starting operator at length $n")
    flush(stdout)
    if get(data, :construction, nothing) == "quartic-taylor-sos"
        read_factor = row -> begin
            witness = joinpath(dirname(source), String(row.witness_file))
            text = read(witness, String)
            @assert bytes2hex(sha256(codeunits(text))) == row.witness_sha256
            return JSON3.read(text)
        end
        seed_progress = index -> begin
            println("Checked starting quartic factor $index/3")
            flush(stdout)
        end
        checked_seed = PM.verify_quartic_seed_data(data, read_factor; progress = seed_progress)
        @assert checked_seed.candidate == B && checked_seed.error == epsilon
    elseif get(data, :exact_compression, false)
        @assert PM.verify_exact_seed(PM.normalized_polynomial(seed), B, data)
        @assert epsilon == 0
    else
        initial = PM.initial_error(PM.normalized_polynomial(seed), B)
        @assert epsilon >= initial.error
    end
    RR = Nemo.ArbField(192)
    logscale = RR(0)
    reference_norm = PM.initial_operator_upper(seed)
    normalizations, local_errors, local_bounds = MN.F[], MN.F[], MN.F[]
    norm_upper, beta = MN.F(0), MN.F(0)
    physical = sqrt.([Float64(PM.WEIGHTS[t]/PM.WEIGHTS[s]) for s in 1:5, t in 1:5])
    result = Dict{String,Any}(
        "status"=>"certified-polynomial-squaring-prefix",
        "requested_steps_complete"=>false,
        "target_beta_reached"=>false,
        "target_thermal_trace_certified"=>false,
        "uniform_gap_certified"=>false,
        "source_file"=>source,
        "source_sha256"=>bytes2hex(sha256(read(source))),
        "shift"=>string(shift),
        "gram_steps"=>gram,
        "reflected_outputs"=>reflect_output,
        "conserve_charge"=>conserve_charge,
        "error_rounding_bits"=>error_round_bits,
        "analytic_error_propagation"=>true,
        "runs"=>Any[],
    )
    first_step = 1
    if !isnothing(resume_from)
        # A generator may reuse saved candidates, but the final verifier still
        # replays every old and new local witness from the original seed.
        previous = JSON3.read(read(resume_from, String), Dict{String,Any})
        # A saved prefix may be valid even if a longer requested run has not
        # finished. Every retained step must still have its complete witness.
        @assert !isempty(previous["runs"])
        @assert previous["source_file"] == source &&
                previous["source_sha256"] == result["source_sha256"]
        @assert MN.parse_fraction(previous["shift"]) == shift
        @assert previous["gram_steps"] == gram && previous["reflected_outputs"] == reflect_output
        @assert abspath(dirname(resume_from)) == abspath(dirname(output))
        @assert length(previous["runs"]) < min(steps, power)
        for (j, row) in enumerate(previous["runs"])
            @assert row["step"] == j
            file = joinpath(dirname(resume_from), String(row["witness_file"]))
            @assert bytes2hex(sha256(read(file))) == row["witness_sha256"]
            proof = JSON3.read(read(file, String))
            @assert MN.decode_cores(proof.input) == B
            @assert MN.parse_fraction(proof.previous_error_upper) == epsilon
            kappa = MN.parse_fraction(proof.normalization)
            @assert kappa > 0
            push!(normalizations, kappa)
            push!(
                local_errors,
                MN.parse_fraction(proof.error)+kappa*MN.parse_fraction(
                    proof.canonicalization_error_upper,
                ),
            )
            push!(local_bounds, MN.parse_fraction(proof.operator_norm_upper))
            logscale = 2logscale+log(RR(kappa))
            reference_norm = PM.PG.relative_upper(RR(reference_norm)^2/RR(kappa))
            B = MS.flow_next_candidate(proof)
            epsilon = MN.parse_fraction(proof.next_error_upper)
            norm_upper = MN.parse_fraction(proof.next_hs_norm_upper)
            beta = MN.parse_fraction(proof.current_beta)
            @assert beta == 2seed.delta*big(2)^j
            push!(result["runs"], row)
        end
        first_step = length(previous["runs"])+1
        result["resumed_from"] = resume_from
        result["resumed_from_sha256"] = bytes2hex(sha256(read(resume_from)))
        println(
            "Continuing from $(first_step-1) saved squarings; old witnesses remain in the final replay",
        )
        flush(stdout)
    end
    for j in first_step:min(steps, power)
        started = time()
        println("Starting polynomial squaring $j/$power at length $n")
        flush(stdout)
        gauges = if j == first_step && !isnothing(reuse_square_from)
            nothing
        else
            first_copy = gram ? MN.adjoint_cores(B; weights = PM.WEIGHTS) : B
            numerical = [Float64.(c) .* reshape(physical, 1, 5, 5, 1) for c in first_copy]
            @assert gauge_mode in (:diagonal, :schmidt)
            selected, objective = if gauge_mode == :schmidt
                proposed = GS.schmidt_gauge_proposal(numerical; physical_charges)
                rational = [MS.relative_quantize(g; bits = 64) for g in proposed]
                rational, GS.matrix_gauge_objective(numerical, proposed)
            else
                search = optimize_gauge(numerical; iterations = gauge_iterations)
                MN.rational_gauges(B, search.x), search.final
            end
            println(
                "Gauge search complete: proposed operator/reference norm ratio $(exp(objective)/Float64(reference_norm))",
            )
            flush(stdout)
            selected
        end
        progress =
            row -> begin
                if row.site == 0 || row.site == row.total || row.site % 32 == 0
                    println("Squaring $j: $(row.stage), site $(row.site)/$(row.total)")
                    flush(stdout)
                end
            end
        square = if j == first_step && !isnothing(reuse_square_from)
            input, saved = MS.decode_witness(JSON3.read(read(reuse_square_from, String)))
            @assert input == B && saved.weights == PM.WEIGHTS && saved.adjoint_first == gram
            if saved.error_method == error_method && isnothing(reuse_gauge_power)
                saved
            else
                fixed = if isnothing(reuse_gauge_power)
                    saved
                else
                    @assert 0 <= reuse_gauge_power <= 1
                    @assert first(saved.gauges) isa AbstractMatrix
                    scales = [
                        MS.relative_quantize(
                            vec(sqrt.(sum(abs2, Float64.(g); dims = 2))) .^ (2reuse_gauge_power-1);
                            bits = 64,
                        ) for g in saved.gauges
                    ]
                    MS.rescale_square_gauges(saved, scales)
                end
                MS.square_candidate(
                    B;
                    gauges = fixed.gauges,
                    weights = PM.WEIGHTS,
                    blocksize = isempty(saved.tail.block_bounds) ? 1 : 2,
                    adjoint_first = gram,
                    candidate_arithmetic,
                    error_method,
                    fixed_compression = fixed,
                    progress,
                )
            end
        else
            MS.square_candidate(
                B;
                maxdim = output_dimension,
                bits = compression_bits,
                gauges,
                weights = PM.WEIGHTS,
                blocksize = tail_blocksize,
                adjoint_first = gram,
                candidate_arithmetic,
                error_method,
                project_bridge,
                physical_charges,
                progress,
            )
        end
        gauge_refinement = nothing
        if !isnothing(gauge_refiner)
            gauge_refinement = gauge_refiner(B, square, reference_norm, epsilon)
            @assert gauge_refinement.certificate.cores == square.cores
            square = gauge_refinement.certificate
        end
        proofpath = splitext(output)[1]*"-step$(j)-witness.json"
        if j == first_step && !isnothing(reuse_square_from)
            @assert abspath(proofpath) != abspath(reuse_square_from)
        end
        proof = MS.witness_encoding(B, square)
        # Preserve the expensive local witness before normalization or any
        # later contraction. The result file is finalized only after checks.
        open(proofpath, "w") do io
            JSON3.write(io, proof)
        end
        println("Replaying local square certificate")
        flush(stdout)
        checked = MS.verify_square_certificate(
            B,
            square;
            details = true,
            residual_arithmetic = :interval,
            progress,
        )
        b = MN.relative_rational_upper(min(checked.operator_norm_upper, checked.hs_norm_upper))
        C = square.cores
        # The stored left-map bounds control the complete HS norm, avoiding
        # interval wrapping in a long dense double-layer contraction.
        kappa = MN.relative_rational_upper(checked.output_hs_norm_upper)
        @assert kappa > 0
        if reflect_output
            next = MN.reflect_cores(C)
            next[1] = next[1]/kappa
            recanonical = (; error = MN.F(0))
            norm_upper = MN.F(1)
        else
            scaled = copy(C)
            scaled[1] = scaled[1]/kappa
            next = candidate_mpo(scaled; maxdim = retained_dimension, bits = 48)
            recanonical = PM.initial_error(scaled, next)
            norm_upper = PM.PG.relative_upper(sqrt(recanonical.cc))
        end
        coefficient = MS.gram_error_coefficient(epsilon, b, reference_norm)
        next_epsilon = (square.error+epsilon*coefficient)/kappa+recanonical.error
        if !isnothing(error_round_bits)
            next_epsilon = MN.relative_rational_upper(next_epsilon; bits = error_round_bits)
        end
        push!(normalizations, kappa)
        push!(local_errors, square.error+kappa*recanonical.error)
        push!(local_bounds, b)
        logscale = 2logscale+log(RR(kappa))
        beta = 2seed.delta*big(2)^j
        logtrace =
            -RR(shift-1)*RR(beta)*n +
            n*big(2)^j*log(RR(5)) +
            2logscale +
            2log(RR(norm_upper+next_epsilon))
        trace_upper = PM.PG.rational_upper(logtrace)
        proof["normalization"] = string(kappa)
        proof["reflected_output"] = reflect_output
        if reflect_output
            # Reflection and the exact scalar already determine every entry.
            # Avoid serializing a second copy of all output tensor coefficients.
            proof["next_candidate_encoding"] = "reflected-output"
        else
            proof["next_candidate"] = MN.encode_cores(next)
        end
        proof["canonicalization_error_upper"] = string(recanonical.error)
        proof["previous_error_upper"] = string(epsilon)
        proof["next_error_upper"] = string(next_epsilon)
        proof["operator_norm_upper"] = string(b)
        proof["next_hs_norm_upper"] = string(norm_upper)
        proof["normalization_norm_method"] = "left-map-product"
        proof["current_beta"] = string(beta)
        proof["shifted_log_trace_upper"] = string(trace_upper)
        open(proofpath, "w") do io
            JSON3.write(io, proof)
        end
        row = Dict{String,Any}(
            "step"=>j,
            "current_beta"=>string(beta),
            "current_thermal_trace_upper_certified"=>true,
            "propagated_error_upper"=>string(next_epsilon),
            "propagated_error_display"=>Float64(next_epsilon),
            "relative_squaring_error_display"=>Float64(square.error/kappa),
            "canonicalization_error_display"=>Float64(recanonical.error),
            "squaring_error_amplification_display"=>Float64(coefficient/kappa),
            "shifted_log_trace_upper"=>string(trace_upper),
            "shifted_log_trace_upper_display"=>Float64(trace_upper),
            "witness_file"=>basename(proofpath),
            "witness_sha256"=>bytes2hex(sha256(read(proofpath))),
            "seconds"=>time()-started,
        )
        if !isnothing(gauge_refinement)
            row["selected_gauge_power"] = gauge_refinement.power
            row["numerical_gauge_screen_certified"] = false
            row["numerical_gauge_screen"] = gauge_refinement.screen
        end
        if n == 4
            spectrum = [
                (div(s*(s+1)-a*(a+1)-c*(c+1), 2), 2s+1) for a in 0:4 for c in 0:4 for
                s in abs(a-c):(a+c)
            ]
            exact_log = log(sum(deg*exp(-RR(beta)*e) for (e, deg) in spectrum))-RR(shift)*RR(beta)*n
            @assert RR(trace_upper) > exact_log
            row["upper_over_exact_trace_display"] = Float64(exp(RR(trace_upper)-exact_log))
        end
        push!(result["runs"], row)
        write_result(output, result)
        println(
            "step=$j beta=$(Float64(beta)) error=$(Float64(next_epsilon)) amplification=$(Float64(coefficient/kappa)) logZupper=$(Float64(trace_upper)) seconds=$(time()-started)",
        )
        flush(stdout)
        B, epsilon = next, next_epsilon
        reference_norm = PM.PG.relative_upper(RR(reference_norm)^2/RR(kappa))
    end
    result["requested_steps_complete"] = true
    if gram && !isempty(normalizations)
        bounds = MS.backward_operator_bounds(normalizations, local_errors, norm_upper)
        for j in eachindex(local_bounds)
            bounds[j] = min(bounds[j], local_bounds[j])
        end
        errors = MS.propagate_gram_errors(
            initial_epsilon,
            normalizations,
            local_errors,
            bounds;
            round_bits = error_round_bits,
            reference_bounds = PM.reference_operator_bounds(seed, normalizations),
        )
        refined_log =
            -RR(shift-1)*RR(beta)*n +
            n*big(2)^length(normalizations)*log(RR(5)) +
            2logscale +
            2log(RR(norm_upper+last(errors)))
        result["backward_operator_bounds"] = string.(bounds)
        result["backward_propagated_errors"] = string.(errors)
        result["backward_final_error_display"] = Float64(last(errors))
        result["backward_shifted_log_trace_upper"] = string(PM.PG.rational_upper(refined_log))
        println("Backward norm refinement: final error <= $(Float64(last(errors)))")
        flush(stdout)
    end
    result["target_beta_reached"] = steps >= power
    result["target_thermal_trace_certified"] = steps >= power
    write_result(output, result)
end

if abspath(PROGRAM_FILE) == @__FILE__
    @assert length(ARGS) == 2
    flow_probe(ARGS...)
end
