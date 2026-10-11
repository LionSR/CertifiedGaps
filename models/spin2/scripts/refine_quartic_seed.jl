include(joinpath(@__DIR__, "gauge_error_screen.jl"))

"Refine a quartic certificate by rescaling gauges, keeping its candidate exact."
function refine_quartic_seed(source, output; powers = 0.0:0.05:0.7)
    @assert abspath(source) != abspath(output)
    @assert abspath(dirname(source)) == abspath(dirname(output))
    data = JSON3.read(read(source, String), Dict{String,Any})
    @assert data["construction"] == "quartic-taylor-sos" && data["degree"] == 4
    seed = PM.PG.taylor_seed(data["length"], MN.parse_fraction(data["beta"]), data["power"], 4)
    @assert seed.remainder == MN.parse_fraction(data["remainder"])
    @assert length(data["factors"]) == 3
    factors = PM.quartic_sos_factors(seed).factors
    sigma = MN.F(1, big(5)^div(seed.n, 2))
    outputs, errors, rows = Vector{Array{MN.F,4}}[], MN.F[], Any[]
    for (index, old) in enumerate(data["factors"])
        @assert old["index"] == index
        witness = joinpath(dirname(source), old["witness_file"])
        @assert bytes2hex(sha256(read(witness))) == old["witness_sha256"]
        raw = JSON3.read(read(witness, String))
        input, saved = MS.decode_witness(raw)
        exact_input = get(raw, :exact_input, nothing)
        initial_epsilon = if isnothing(exact_input)
            PM.initial_error(PM.normalized_polynomial(factors[index]), input).error
        else
            @assert PM.verify_exact_seed(PM.normalized_polynomial(factors[index]), input, exact_input)
            MN.F(0)
        end
        println("Screening quartic factor $index/3")
        flush(stdout)
        refined = refine_square_gauge(
            input,
            saved,
            PM.initial_operator_upper(factors[index]),
            initial_epsilon;
            powers,
        )
        certificate, power, screen = refined.certificate, refined.power, refined.screen
        path = splitext(output)[1]*"-factor$(index)-witness.json"
        witness_data = MS.witness_encoding(input, certificate)
        isnothing(exact_input) || (witness_data["exact_input"] = exact_input)
        open(path, "w") do io
            JSON3.write(io, witness_data)
        end
        checked = PM.verify_quartic_factor(seed, index, input, certificate; exact_input)
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
                "selected_gauge_power"=>power,
                "numerical_screen_certified"=>false,
                "numerical_screen_normalization"=>"unscaled Gram error",
                "numerical_screen"=>screen,
            ),
        )
        println(
            "Refined factor $index: power=$power normalized error<=$(Float64(checked.error/sigma))",
        )
        flush(stdout)
    end
    candidate = MN.decode_cores(data["candidate"])
    checked = PM.quartic_sos_seed_error(seed, outputs, errors, candidate)
    data["factors"] = rows
    data["error_upper"] = string(checked.error)
    data["factor_error_upper"] = string(checked.factor_error)
    data["final_compression_error_upper"] = string(checked.compression_error)
    data["refined_from"] = source
    data["refined_from_sha256"] = bytes2hex(sha256(read(source)))
    data["refinement_runner_sha256"] = bytes2hex(sha256(read(@__FILE__)))
    data["unchanged_candidate"] = true
    open(output, "w") do io
        JSON3.write(io, data)
    end
    println("Refined quartic seed: error<=$(Float64(checked.error))")
    return data
end

if abspath(PROGRAM_FILE) == @__FILE__
    @assert length(ARGS) == 2
    refine_quartic_seed(ARGS...)
end
