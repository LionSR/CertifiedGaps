include(joinpath(@__DIR__, "polynomial_flow_probe.jl"))
using ThermalMPOGap.Results: write_result

"Numerical screening of diagonal rescalings of a saved full-matrix gauge."
function square_gauge_error_screen(
    B,
    proof;
    reference,
    initial_error,
    kappa = 1,
    powers = 0.15:0.05:0.6,
    precision = 192,
    independent_operator_upper = nothing,
)
    @assert proof.adjoint_first && proof.weights == PM.WEIGHTS
    @assert first(proof.gauges) isa AbstractMatrix
    n = length(B)
    RR = Nemo.ArbField(precision)
    gauged = MN.gauge_transform(
        MN.adjoint_cores(B; weights = PM.WEIGHTS),
        proof.gauges;
        interval_precision = precision,
    )
    physical = sqrt.([Float64(PM.WEIGHTS[t]/PM.WEIGHTS[s]) for s in 1:5, t in 1:5])
    numeric = [Float64.(c) .* reshape(physical, 1, 5, 5, 1) for c in gauged]
    residual_buckets, projection_buckets = Vector{Float64}[], Vector{Float64}[]
    bridge = ones(MN.F, 1, 1, 1)
    for j in 1:(n-1)
        Z = MS.square_bridge(bridge, gauged[j], B[j]; field = RR)
        chi, d, _, r, _ = size(Z)
        M = reshape(Z, chi*d*d, r*r)
        k = size(proof.cores[j], 4)
        U, R = reshape(proof.cores[j], chi*d*d, k), reshape(proof.bridges[j], k, r*r)
        E = M-MS.matrix_product(U, R, RR)
        metric = MN.hs_row_metric(chi, d, PM.WEIGHTS)
        # Midpoints are used only to choose a prospective gauge. No numerical
        # bucket or screening result is accepted as a certificate input.
        weighted = Float64.(E) .* sqrt.(Float64.(metric))
        projected = Float64.(MS.matrix_product(transpose(U .* metric), E, RR))
        push!(residual_buckets, vec(sum(abs2, reshape(weighted, chi*d*d, r, r); dims = (1, 3))))
        push!(projection_buckets, vec(sum(abs2, reshape(projected, k, r, r); dims = (1, 3))))
        bridge = proof.bridges[j]
        if j % 32 == 0
            println("Numerical gauge-error screen: residual $j/$(n-1)")
            flush(stdout)
        end
    end
    base_gauges = [Float64.(g) for g in proof.gauges]
    row_scales = [vec(sqrt.(sum(abs2, g; dims = 2))) for g in base_gauges]
    local_hs = Float64.(proof.tail.local_hs_bounds)
    alpha = [Float64(row.local_prefix_bound) for row in proof.rows]
    hs_tail = ones(n+1)
    for j in n:-1:1
        hs_tail[j] = local_hs[j]*hs_tail[j+1]
    end
    reference, kappa, initial_error = Float64(reference), Float64(kappa), Float64(initial_error)
    @assert reference > 0 && kappa > 0 && initial_error >= 0
    rows = Any[]
    for power in powers
        scales = [s .^ (2power-1) for s in row_scales]
        op_tail = ones(n+1)
        for j in n:-1:1
            C = numeric[j] .* reshape(scales[j], :, 1, 1, 1) ./ reshape(scales[j+1], 1, 1, 1, :)
            local_op = opnorm(MN.local_operator(C))
            cap = hs_tail[j]*norm(base_gauges[j] .* scales[j])
            op_tail[j] = min(local_op*op_tail[j+1], cap)
        end
        error, triangle = 0.0, 0.0
        for j in (n-1):-1:1
            delta = sqrt(sum(residual_buckets[j] ./ scales[j+1] .^ 2))
            theta = min(sqrt(sum(projection_buckets[j] ./ scales[j+1] .^ 2)), alpha[j]*delta)
            tail = op_tail[j+1]*hs_tail[j+1]
            error = sqrt((delta*tail)^2+(alpha[j]*error)^2+2theta*tail*error)
            triangle = delta*tail+alpha[j]*triangle
        end
        bound = min(op_tail[1], hs_tail[1])
        if !isnothing(independent_operator_upper)
            bound = min(bound, Float64(independent_operator_upper))
        end
        coefficient = min(2reference+initial_error, 2bound+initial_error, reference+bound)
        row = Dict(
            "power"=>power,
            "operator_over_reference_estimate"=>bound/reference,
            "projection_relative_error_estimate"=>error/kappa,
            "triangle_relative_error_estimate"=>triangle/kappa,
            "propagated_error_estimate"=>(error+coefficient*initial_error)/kappa,
        )
        push!(rows, row)
        println(
            "power=$power projected error estimate=$(row["propagated_error_estimate"]), operator/reference=$(bound/reference)",
        )
        flush(stdout)
    end
    return rows
end

"Propose a gauge rescaling numerically, then certify it with frozen output cores."
function refine_square_gauge(input, saved, reference, initial_error; powers = 0.0:0.05:0.7)
    independent = (; gauges = saved.gauges, bounds = saved.tail.local_operator_bounds)
    if hasproperty(saved, :operator_certificate) &&
       !isnothing(saved.operator_certificate) &&
       prod(saved.operator_certificate.bounds) < prod(independent.bounds)
        independent = saved.operator_certificate
    end
    screen = square_gauge_error_screen(
        input,
        saved;
        reference,
        initial_error,
        powers,
        independent_operator_upper = prod(independent.bounds),
    )
    best = argmin(row -> row["propagated_error_estimate"], screen)
    power = best["power"]
    scales = [
        MS.relative_quantize(vec(sqrt.(sum(abs2, Float64.(g); dims = 2))) .^ (2power-1); bits = 64) for g in saved.gauges
    ]
    fixed = MS.rescale_square_gauges(saved, scales)
    certificate = MS.square_candidate(
        input;
        gauges = fixed.gauges,
        weights = PM.WEIGHTS,
        blocksize = 1,
        adjoint_first = true,
        candidate_arithmetic = :interval,
        error_method = :projection,
        fixed_compression = fixed,
    )
    @assert certificate.cores == saved.cores
    certificate = merge(certificate, (; operator_certificate = independent))
    return (; certificate, power, screen)
end

function gauge_error_screen(source, output; powers = 0.15:0.05:0.6, precision = 192)
    source_text = read(source, String)
    source_hash = bytes2hex(sha256(codeunits(source_text)))
    flow = JSON3.read(source_text)
    @assert !isempty(flow.runs)
    step = length(flow.runs)
    @assert last(flow.runs).step == step
    witness = joinpath(dirname(source), String(last(flow.runs).witness_file))
    witness_text = read(witness, String)
    witness_hash = bytes2hex(sha256(codeunits(witness_text)))
    @assert witness_hash == last(flow.runs).witness_sha256
    data = JSON3.read(witness_text)
    B, proof = MS.decode_witness(data)
    initial = JSON3.read(read(String(flow.source_file), String))
    @assert Int(initial.length) == length(B)
    seed = PM.scalar_seed(initial)
    kappas = MN.F[]
    for row in flow.runs[1:(end-1)]
        prior_path = joinpath(dirname(source), String(row.witness_file))
        prior_text = read(prior_path, String)
        @assert bytes2hex(sha256(codeunits(prior_text))) == row.witness_sha256
        prior = JSON3.read(prior_text)
        push!(kappas, MN.parse_fraction(prior.normalization))
    end
    rows = square_gauge_error_screen(
        B,
        proof;
        reference = last(PM.reference_operator_bounds(seed, kappas)),
        initial_error = MN.parse_fraction(data.previous_error_upper),
        kappa = MN.parse_fraction(data.normalization),
        powers,
        precision,
    )
    result = Dict(
        "status"=>"numerical-gauge-error-screen-only",
        "certified"=>false,
        "uniform_gap_certified"=>false,
        "source_file"=>source,
        "source_sha256"=>source_hash,
        "witness_file"=>witness,
        "witness_sha256"=>witness_hash,
        "step"=>step,
        "rows"=>rows,
    )
    write_result(output, result)
    return result
end

if abspath(PROGRAM_FILE) == @__FILE__
    @assert length(ARGS) == 2
    gauge_error_screen(ARGS...)
end
