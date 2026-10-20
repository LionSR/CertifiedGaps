include(joinpath(@__DIR__, "polynomial_flow_probe.jl"))

"Build a fourth-order seed using three existing quadratic Gram certificates."
function quartic_seed_probe(
    output;
    n = 4,
    beta = 4,
    power = 8,
    maxdim = 32,
    input_dimension = min(maxdim, 34),
    bits = 64,
    exact_inputs = false,
    condition_exact_inputs = false,
)
    @assert maxdim >= 1 && input_dimension >= 1
    @assert !condition_exact_inputs || exact_inputs
    seed = PM.PG.taylor_seed(n, beta, power, 4)
    physical = sqrt.([Float64(PM.WEIGHTS[t]/PM.WEIGHTS[s]) for s in 1:5, t in 1:5])
    outputs, errors, rows = Vector{Array{MN.F,4}}[], MN.F[], Any[]
    for (index, factor) in enumerate(PM.quartic_sos_factors(seed).factors)
        println("Quartic seed: factor $index/3 at N=$n")
        flush(stdout)
        original = PM.normalized_polynomial(factor)
        input, exact_input = if exact_inputs
            right = PM.exact_compress_right(original)
            left = PM.exact_compress_right(MN.reflect_cores(right.cores))
            exact = MN.reflect_cores(left.cores)
            proof = Dict{String,Any}(
                "right_output"=>MN.encode_cores(right.cores),
                "right_transforms"=>PM.encode_transforms(right.transforms),
                "left_transforms"=>PM.encode_transforms(left.transforms),
            )
            if condition_exact_inputs
                conditioned = PM.condition_right_exact(exact; precision = 384, bits = 80)
                proof["exact_conditioning"] = true
                proof["unconditioned_candidate"] = MN.encode_cores(exact)
                proof["conditioning_transforms"] = PM.encode_transforms(conditioned.transforms)
                exact = conditioned.cores
            end
            @assert PM.verify_exact_seed(original, exact, proof)
            exact, proof
        else
            candidate_mpo(original; maxdim = input_dimension, bits), nothing
        end
        numerical = [
            Float64.(c) .* reshape(physical, 1, 5, 5, 1) for
            c in MN.adjoint_cores(input; weights = PM.WEIGHTS)
        ]
        gauges = [MS.relative_quantize(g; bits = 64) for g in GS.schmidt_gauge_proposal(numerical)]
        progress =
            row -> begin
                if row.site == row.total || row.site % 32 == 0
                    println("Quartic factor $index: $(row.stage), site $(row.site)/$(row.total)")
                    flush(stdout)
                end
            end
        certificate = MS.square_candidate(
            input;
            maxdim,
            bits,
            gauges,
            weights = PM.WEIGHTS,
            blocksize = 1,
            adjoint_first = true,
            candidate_arithmetic = :interval,
            error_method = :projection,
            project_bridge = true,
            progress,
        )
        path = splitext(output)[1]*"-factor$(index)-witness.json"
        witness_data = MS.witness_encoding(input, certificate)
        isnothing(exact_input) || (witness_data["exact_input"] = exact_input)
        open(path, "w") do io
            JSON3.write(io, witness_data)
        end
        checked = PM.verify_quartic_factor(seed, index, input, certificate; progress, exact_input)
        push!(outputs, certificate.cores)
        push!(errors, checked.error)
        push!(
            rows,
            Dict(
                "index"=>index,
                "witness_file"=>basename(path),
                "witness_sha256"=>bytes2hex(sha256(read(path))),
                "factor_initial_error_upper"=>string(checked.initial_error),
                "gram_error_upper"=>string(checked.error),
                "input_dimension"=>maximum(size(c, 4) for c in input),
            ),
        )
    end
    println("Quartic seed: certifying final compression")
    flush(stdout)
    combined = PM.quartic_sos_assembly(seed, outputs)
    candidate = candidate_mpo(combined; maxdim, bits)
    checked = PM.quartic_sos_seed_error(seed, outputs, errors, candidate)
    result = Dict(
        "construction"=>"quartic-taylor-sos",
        "length"=>n,
        "beta"=>string(seed.beta),
        "power"=>power,
        "degree"=>4,
        "remainder"=>string(seed.remainder),
        "candidate"=>MN.encode_cores(candidate),
        "error_upper"=>string(checked.error),
        "factor_error_upper"=>string(checked.factor_error),
        "final_compression_error_upper"=>string(checked.compression_error),
        "precision_bits"=>192,
        "maxdim"=>maxdim,
        "input_dimension"=>maximum(row["input_dimension"] for row in rows),
        "exact_inputs"=>exact_inputs,
        "condition_exact_inputs"=>condition_exact_inputs,
        "factors"=>rows,
        "initial_error_certified"=>true,
        "thermal_trace_certified"=>false,
        "uniform_gap_certified"=>false,
        "runner_sha256"=>bytes2hex(sha256(read(@__FILE__))),
    )
    open(output, "w") do io
        JSON3.write(io, result)
    end
    println("Quartic seed complete: N=$n error<=$(Float64(checked.error))")
    flush(stdout)
    return result
end

if abspath(PROGRAM_FILE) == @__FILE__
    @assert length(ARGS) == 1
    quartic_seed_probe(only(ARGS))
end
