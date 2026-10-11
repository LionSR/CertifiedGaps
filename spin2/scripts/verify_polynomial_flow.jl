include(joinpath(@__DIR__, "verify_polynomial_seed.jl"))
using SHA
import Nemo
const MS = PM.MS

function verify_flow(path; residual_arithmetic = :interval)
    result = JSON3.read(read(path, String))
    source = String(result.source_file)
    @assert bytes2hex(sha256(read(source))) == result.source_sha256
    initial_epsilon = verify_seed(source; residual_arithmetic)
    data = JSON3.read(read(source, String))
    n, power = Int(data.length), Int(data.power)
    seed = PM.scalar_seed(data)
    B = MN.decode_cores(data.candidate)
    epsilon = initial_epsilon
    RR = Nemo.ArbField(192)
    logscale = RR(0)
    reference_norm = PM.initial_operator_upper(seed)
    shift = MN.parse_fraction(result.shift)
    kappas, deltas, local_bounds = MN.F[], MN.F[], MN.F[]
    norm_upper, beta = MN.F(0), MN.F(0)
    all_gram = true
    for (j, row) in enumerate(result.runs)
        @assert row.step == j && j <= power
        witness_path = joinpath(dirname(path), String(row.witness_file))
        @assert bytes2hex(sha256(read(witness_path))) == row.witness_sha256
        proof = JSON3.read(read(witness_path, String))
        input, certificate = MS.decode_witness(proof)
        @assert input == B && certificate.weights == PM.WEIGHTS
        progress =
            row -> begin
                if row.site % 32 == 0 || row.site == row.total
                    println("Replay step $j: $(row.stage), site $(row.site)/$(row.total)")
                    flush(stdout)
                end
            end
        checked_local = MS.verify_square_certificate(
            input,
            certificate;
            details = true,
            residual_arithmetic,
            progress,
        )
        all_gram &= certificate.adjoint_first
        b = MN.parse_fraction(proof.operator_norm_upper)
        @assert b >= min(checked_local.operator_norm_upper, checked_local.hs_norm_upper)
        kappa = MN.parse_fraction(proof.normalization)
        @assert kappa > 0
        next = MS.flow_next_candidate(proof; output_cores = certificate.cores)
        rho = MN.parse_fraction(proof.canonicalization_error_upper)
        @assert rho >= 0
        if get(proof, :reflected_output, false)
            expected = MN.reflect_cores(certificate.cores)
            expected[1] = expected[1]/kappa
            @assert next == expected
            cc =
                get(proof, :normalization_norm_method, "contraction") == "left-map-product" ?
                nothing : MS.hs_inner_ball(next, next; weights = PM.WEIGHTS)
        else
            scaled = copy(certificate.cores)
            scaled[1] = scaled[1]/kappa
            canonical = PM.initial_error(scaled, next)
            @assert rho >= canonical.error
            cc = canonical.cc
        end
        @assert MN.parse_fraction(proof.previous_error_upper) == epsilon
        next_epsilon = MN.parse_fraction(proof.next_error_upper)
        coefficient = MS.gram_error_coefficient(epsilon, b, reference_norm)
        @assert next_epsilon >= (certificate.error+epsilon*coefficient)/kappa+rho
        norm_upper = MN.parse_fraction(proof.next_hs_norm_upper)
        @assert norm_upper > 0
        if isnothing(cc)
            @assert norm_upper*kappa >= checked_local.output_hs_norm_upper
        else
            @assert cc < RR(norm_upper)^2
        end
        push!(kappas, kappa)
        push!(deltas, certificate.error+kappa*rho)
        push!(local_bounds, b)
        logscale = 2logscale+log(RR(kappa))
        beta = 2seed.delta*big(2)^j
        @assert MN.parse_fraction(proof.current_beta) == beta
        logtrace =
            -RR(shift-1)*RR(beta)*n+n*big(2)^j*log(RR(5))+2logscale+2log(
                RR(norm_upper+next_epsilon),
            )
        bound = MN.parse_fraction(proof.shifted_log_trace_upper)
        @assert RR(bound) > logtrace && MN.parse_fraction(row.shifted_log_trace_upper) == bound
        B, epsilon = next, next_epsilon
        reference_norm = PM.PG.relative_upper(RR(reference_norm)^2/RR(kappa))
    end
    if haskey(result, :backward_operator_bounds)
        @assert all_gram
        bounds = MS.backward_operator_bounds(kappas, deltas, norm_upper)
        for j in eachindex(local_bounds)
            bounds[j] = min(bounds[j], local_bounds[j])
        end
        @assert MN.parse_fraction.(result.backward_operator_bounds) == bounds
        errors = MS.propagate_gram_errors(
            initial_epsilon,
            kappas,
            deltas,
            bounds;
            round_bits = get(result, :error_rounding_bits, nothing),
            reference_bounds = get(result, :analytic_error_propagation, false) ?
                               PM.reference_operator_bounds(seed, kappas) : nothing,
        )
        @assert MN.parse_fraction.(result.backward_propagated_errors) == errors
        logtrace =
            -RR(shift-1)*RR(beta)*n+n*big(2)^length(kappas)*log(RR(5))+2logscale+2log(
                RR(norm_upper+last(errors)),
            )
        @assert RR(MN.parse_fraction(result.backward_shifted_log_trace_upper)) > logtrace
    end
    println(
        "Verified $(length(kappas)) thermal squarings at N=$n; final beta=$(Float64(beta)); target beta reached=$(length(kappas)==power)",
    )
    return (; steps = length(kappas), beta, target_reached = length(kappas)==power)
end

if abspath(PROGRAM_FILE) == @__FILE__
    for path in ARGS
        verify_flow(path)
    end
end
