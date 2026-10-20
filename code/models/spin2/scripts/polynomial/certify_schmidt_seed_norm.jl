include(joinpath(@__DIR__, "polynomial_flow_probe.jl"))

function certify_schmidt_seed_norm(source, output; bits = 64, precision = 192)
    started = time()
    data = JSON3.read(read(source, String))
    seed = PM.scalar_seed(data)
    @assert seed.remainder == MN.parse_fraction(data.remainder)
    B = MN.decode_cores(data.candidate)
    cores = MN.adjoint_cores(B; weights = PM.WEIGHTS)
    physical = sqrt.([Float64(PM.WEIGHTS[t]/PM.WEIGHTS[s]) for s in 1:5, t in 1:5])
    numerical = [Float64.(c) .* reshape(physical, 1, 5, 5, 1) for c in cores]
    proposed = GS.schmidt_gauge_proposal(numerical)
    gauges = [MS.relative_quantize(g; bits) for g in proposed]
    certificate = MN.local_product_witness(cores, gauges; precision, weights = PM.WEIGHTS)
    reference = PM.initial_operator_upper(seed)
    result = Dict(
        "status"=>"exact-rational-candidate-MPO-norm-witness",
        "source_file"=>source,
        "source_sha256"=>bytes2hex(sha256(read(source))),
        "source_adjoint"=>true,
        "length"=>length(cores),
        "cores"=>MN.encode_cores(cores),
        "gauges"=>MN.encode_gauges(gauges),
        "physical_weights"=>string.(PM.WEIGHTS),
        "precision_bits"=>precision,
        "local_bounds"=>string.(certificate.bounds),
        "upper"=>string(certificate.upper),
        "analytic_reference_upper"=>string(reference),
        "operator_over_reference_upper_display"=>Float64(certificate.upper/reference),
        "seconds"=>time()-started,
        "thermal_trace_certified"=>false,
        "uniform_gap_certified"=>false,
    )
    open(output, "w") do io
        JSON3.write(io, result)
    end
    println(
        "Certified N=$(length(cores)) candidate operator/reference upper ratio $(Float64(certificate.upper/reference)); seconds=$(time()-started)",
    )
    return result
end

if abspath(PROGRAM_FILE) == @__FILE__
    @assert length(ARGS) == 2
    certify_schmidt_seed_norm(ARGS...)
end
