include(joinpath(@__DIR__, "../polynomial/verify_polynomial_flow.jl"))
import CertifiedGaps
using CertifiedGaps.Results: write_result
const MP = CertifiedGaps.CertifiedMPO.Products
const AT = CertifiedGaps.CertifiedMPO.AdjointTrace

"Load a fully replayed reflected Gram flow and sharpen its primal error bounds."
function adjoint_primal(path)
    verified = verify_flow(path)
    data = JSON3.read(read(path, String))
    initial = JSON3.read(read(String(data.source_file), String))
    seed = PM.scalar_seed(initial)
    B = [MN.decode_cores(initial.candidate)]
    gauges = Any[]
    kappas, deltas, local_bounds = MN.F[], MN.F[], MN.F[]
    for row in data.runs
        file = joinpath(dirname(path), String(row.witness_file))
        proof = JSON3.read(read(file, String))
        @assert proof.adjoint_first && proof.reflected_output
        @assert MN.parse_fraction(proof.canonicalization_error_upper) == 0
        push!(B, MS.flow_next_candidate(proof))
        push!(kappas, MN.parse_fraction(proof.normalization))
        push!(deltas, MN.parse_fraction(proof.error))
        push!(local_bounds, MN.parse_fraction(proof.operator_norm_upper))
        push!(gauges, haskey(proof, :gauges) ? MN.decode_gauges(proof.gauges) : nothing)
    end
    @assert !isempty(kappas)
    reference = PM.reference_operator_bounds(seed, kappas)
    coefficients = vcat(local_bounds, reference[end])
    errors = MS.propagate_gram_errors(
        MN.parse_fraction(initial.error_upper),
        kappas,
        deltas,
        coefficients,
        round_bits = get(data, :error_rounding_bits, nothing),
        reference_bounds = reference,
    )
    candidate_operator_bounds =
        vcat(min.(local_bounds, reference[1:(end-1)]+errors[1:(end-1)]), reference[end]+errors[end])
    return (;
        B,
        kappas,
        errors,
        gauges,
        candidate_operator_bounds,
        seed,
        shift = MN.parse_fraction(data.shift),
        verified,
    )
end

"Reflection turns a left-to-right product compression into a right-to-left dual."
function dual_product_inputs(B, G)
    first = MN.reflect_cores(B)
    second = PM.linear_combination([G, MN.adjoint_cores(G; weights = PM.WEIGHTS)], MN.F[1, 1])
    return first, second
end

function dual_product_output(proof, kappa)
    G = MN.reflect_cores(proof.cores)
    G[1] = G[1]/kappa
    return G
end

"Dyadic gauges keep odd denominators on outgoing indices, outside bridge sums."
function reflected_dual_gauges(gauges; bits = 40)
    isnothing(gauges) && return nothing
    if first(gauges) isa AbstractMatrix
        return [MS.relative_quantize(transpose(inv(Float64.(g))); bits) for g in reverse(gauges)]
    end
    result = [MS.relative_quantize(Float64.(inv.(g)); bits) for g in reverse(gauges)]
    @assert first(result) == last(result) == MN.F[1]
    @assert all(g -> all(>(0), g), result)
    return result
end

"Contract the signed center; no full physical operator or four-layer MPO is formed."
function adjoint_center(primal, G; precision = 192)
    (; B, kappas, seed) = primal
    weights = PM.WEIGHTS
    cc = MS.hs_inner_ball(last(B), last(B); weights, precision)
    center =
        cc + 2(
            MS.hs_inner_ball(first(G), PM.normalized_polynomial(seed); weights, precision) -
            MS.hs_inner_ball(first(G), first(B); weights, precision)
        )
    terms = String[]
    for j in eachindex(kappas)
        reflected = MN.reflect_cores(G[j+1])
        overlap =
            MP.product_inner_ball(
                reflected,
                MN.adjoint_cores(B[j]; weights),
                B[j];
                weights,
                precision,
            )/parent(cc)(kappas[j])
        residual = overlap - MS.hs_inner_ball(G[j+1], B[j+1]; weights, precision)
        center += 2residual
        push!(terms, string(residual))
        println("Signed residual contraction $j/$(length(kappas)) complete")
        flush(stdout)
    end
    return (; center, cc, terms)
end

function adjoint_evaluation(primal, G, tau, gamma; precision = 192, contracted = nothing)
    (; center, cc, terms) =
        isnothing(contracted) ? adjoint_center(primal, G; precision) : contracted
    lo, hi = -PM.PG.relative_upper(-center), PM.PG.relative_upper(center)
    interval = AT.trace_interval(
        lo,
        hi;
        errors = primal.errors,
        dual_defects = tau,
        dual_operator_bounds = gamma,
        normalizations = primal.kappas,
    )
    @assert interval.upper > 0
    triangle = (PM.PG.relative_upper(sqrt(cc))+last(primal.errors))^2
    RR = parent(cc)
    logscale = RR(0)
    for kappa in primal.kappas
        logscale = 2logscale + log(RR(kappa))
    end
    n, steps = primal.seed.n, length(primal.kappas)
    beta = 2primal.seed.delta*big(2)^steps
    logtrace =
        -RR(primal.shift-1)*RR(beta)*n +
        n*big(2)^steps*log(RR(5)) +
        2logscale +
        log(RR(interval.upper))
    return (;
        interval,
        triangle,
        cc,
        center,
        terms,
        beta,
        logtrace = PM.PG.relative_upper(logtrace),
    )
end

function save_adjoint(path, data)
    write_result(path, data)
end

function mpo_adjoint_flow(
    source,
    output;
    maxdim = 16,
    bits = 64,
    precision = 192,
    reuse_gauges = true,
)
    started = time()
    primal = adjoint_primal(source)
    steps = length(primal.kappas)
    G = Vector{typeof(first(primal.B))}(undef, steps+1)
    G[end] = last(primal.B)
    tau, gamma = zeros(MN.F, steps), zeros(MN.F, steps)
    next_operator_bound = last(primal.candidate_operator_bounds)
    data = Dict{String,Any}(
        "status"=>"multistep-reflected-gram-adjoint-certificate",
        "complete"=>false,
        "source_file"=>source,
        "source_sha256"=>bytes2hex(sha256(read(source))),
        "length"=>primal.seed.n,
        "steps"=>steps,
        "precision_bits"=>precision,
        "max_dual_dimension"=>maxdim,
        "dual_steps"=>Any[],
        "primal_errors"=>string.(primal.errors),
        "source_fully_replayed"=>true,
        "target_thermal_trace_certified"=>false,
        "uniform_gap_certified"=>false,
    )
    save_adjoint(output, data)
    for j in steps:-1:1
        tick = time()
        first, second = dual_product_inputs(primal.B[j], G[j+1])
        # Any positive rational gauge cancels exactly. Round the reciprocal
        # proposal before applying it to avoid large common denominators.
        first_gauges =
            reuse_gauges && !isnothing(primal.gauges[j]) ? reflected_dual_gauges(primal.gauges[j]) :
            nothing
        proof = MP.product_candidate(
            first,
            second;
            maxdim,
            bits,
            precision,
            weights = PM.WEIGHTS,
            first_gauges,
        )
        G[j] = dual_product_output(proof, primal.kappas[j])
        tau[j] = proof.error/primal.kappas[j]
        gg = MS.hs_inner_ball(G[j+1], G[j+1]; weights = PM.WEIGHTS, precision)
        gamma[j] = min(PM.PG.relative_upper(sqrt(gg)), next_operator_bound)
        # The dual defect is also an operator-error bound. This recursion
        # avoids replacing an exponentially small operator norm by an HS norm.
        next_operator_bound = 2primal.candidate_operator_bounds[j]*gamma[j]/primal.kappas[j]+tau[j]
        witness = splitext(output)[1]*"-dual$j-witness.json"
        save_adjoint(witness, MP.encode(first, second, proof))
        push!(
            data["dual_steps"],
            Dict(
                "step"=>j,
                "witness_file"=>basename(witness),
                "witness_sha256"=>bytes2hex(sha256(read(witness))),
                "defect_upper"=>string(tau[j]),
                "operator_norm_upper"=>string(gamma[j]),
                "defect_upper_display"=>Float64(tau[j]),
                "seconds"=>time()-tick,
            ),
        )
        save_adjoint(output, data)
        println(
            "Backward step $j/$steps: defect <= $(Float64(tau[j])), $(round(time()-tick; digits=1)) seconds",
        )
        flush(stdout)
    end
    evaluation = adjoint_evaluation(primal, G, tau, gamma; precision)
    (; interval, triangle, cc, center, terms, beta, logtrace) = evaluation
    merge!(
        data,
        Dict(
            "complete"=>true,
            "current_beta"=>string(beta),
            "target_beta_reached"=>primal.verified.target_reached,
            "polynomial_prefix_trace_certified"=>true,
            "center_ball"=>string(center),
            "signed_residual_balls"=>terms,
            "trace_lower"=>string(interval.lower),
            "trace_upper"=>string(interval.upper),
            "triangle_trace_upper"=>string(triangle),
            "shifted_log_trace_upper"=>string(logtrace),
            "corrected_over_candidate_display"=>Float64(parent(cc)(interval.upper)/cc),
            "triangle_over_candidate_display"=>Float64(parent(cc)(triangle)/cc),
            "cross_remainder_display"=>Float64(interval.cross),
            "quadratic_remainder_display"=>Float64(interval.quadratic),
            "terminal_remainder_display"=>Float64(interval.terminal),
            "seconds"=>time()-started,
        ),
    )
    save_adjoint(output, data)
    println(
        "Full MPO prefix: triangle/candidate=$(data["triangle_over_candidate_display"]), corrected/candidate=$(data["corrected_over_candidate_display"])",
    )
    return data
end

"Replay the primal and dual witnesses and recompute the signed trace enclosure."
function verify_mpo_adjoint_flow(path; refine_output = nothing)
    data = JSON3.read(read(path, String))
    @assert data.complete
    source = String(data.source_file)
    @assert bytes2hex(sha256(read(source))) == data.source_sha256
    primal = adjoint_primal(source)
    steps, precision = length(primal.kappas), Int(data.precision_bits)
    @assert data.steps == steps && data.length == primal.seed.n
    @assert MN.parse_fraction.(data.primal_errors) == primal.errors
    @assert length(data.dual_steps) == steps
    G = Vector{typeof(first(primal.B))}(undef, steps+1)
    G[end] = last(primal.B)
    tau, gamma = zeros(MN.F, steps), zeros(MN.F, steps)
    sharper_gamma = zeros(MN.F, steps)
    next_operator_bound = sharper_next = last(primal.candidate_operator_bounds)
    for (offset, row) in enumerate(data.dual_steps)
        j = steps+1-offset
        @assert row.step == j
        witness = joinpath(dirname(path), String(row.witness_file))
        @assert bytes2hex(sha256(read(witness))) == row.witness_sha256
        first, second, proof = MP.decode(JSON3.read(read(witness, String)))
        expected_first, expected_second = dual_product_inputs(primal.B[j], G[j+1])
        @assert first == expected_first && second == expected_second
        @assert proof.weights == PM.WEIGHTS && MP.verify_product(first, second, proof; precision)
        G[j] = dual_product_output(proof, primal.kappas[j])
        tau[j] = MN.parse_fraction(row.defect_upper)
        @assert tau[j] >= proof.error/primal.kappas[j]
        gamma[j] = MN.parse_fraction(row.operator_norm_upper)
        gg = MS.hs_inner_ball(G[j+1], G[j+1]; weights = PM.WEIGHTS, precision)
        hs_bound = PM.PG.relative_upper(sqrt(gg))
        @assert gamma[j] >= min(hs_bound, next_operator_bound)
        sharper_gamma[j] = min(hs_bound, sharper_next)
        next_operator_bound = 2primal.candidate_operator_bounds[j]*gamma[j]/primal.kappas[j]+tau[j]
        sharper_next = 2primal.candidate_operator_bounds[j]*sharper_gamma[j]/primal.kappas[j]+tau[j]
        println("Verified backward product $j/$steps")
        flush(stdout)
    end
    evaluation = adjoint_evaluation(primal, G, tau, gamma; precision)
    @assert MN.parse_fraction(data.trace_lower) <= evaluation.interval.lower
    @assert MN.parse_fraction(data.trace_upper) >= evaluation.interval.upper
    @assert MN.parse_fraction(data.triangle_trace_upper) >= evaluation.triangle
    @assert MN.parse_fraction(data.shifted_log_trace_upper) >= evaluation.logtrace
    @assert MN.parse_fraction(data.current_beta) == evaluation.beta
    @assert data.target_beta_reached == primal.verified.target_reached
    @assert !data.target_thermal_trace_certified && !data.uniform_gap_certified
    if !isnothing(refine_output)
        refined =
            adjoint_evaluation(primal, G, tau, sharper_gamma; precision, contracted = evaluation)
        result = JSON3.read(JSON3.write(data), Dict{String,Any})
        for row in result["dual_steps"]
            row["operator_norm_upper"] = string(sharper_gamma[row["step"]])
        end
        result["trace_lower"] = string(refined.interval.lower)
        result["trace_upper"] = string(refined.interval.upper)
        result["shifted_log_trace_upper"] = string(refined.logtrace)
        result["corrected_over_candidate_display"] =
            Float64(parent(refined.cc)(refined.interval.upper)/refined.cc)
        result["quadratic_remainder_display"] = Float64(refined.interval.quadratic)
        result["refined_from"] = path
        result["refined_from_sha256"] = bytes2hex(sha256(read(path)))
        result["independently_replayed"] = true
        # Witness paths remain valid if the refinement is saved beside its source.
        @assert abspath(dirname(refine_output)) == abspath(dirname(path))
        save_adjoint(refine_output, result)
        println(
            "Refined dual operator norms: corrected/candidate=$(result["corrected_over_candidate_display"])",
        )
    end
    println(
        "Verified full $(primal.seed.n)-site MPO adjoint prefix at beta=$(Float64(evaluation.beta))",
    )
    return true
end

if abspath(PROGRAM_FILE) == @__FILE__
    if first(ARGS) == "--refine"
        @assert length(ARGS) == 3
        verify_mpo_adjoint_flow(ARGS[2]; refine_output = ARGS[3])
    elseif first(ARGS) == "--verify"
        for path in ARGS[2:end]
            verify_mpo_adjoint_flow(path)
        end
    else
        @assert length(ARGS) in (2, 3)
        mpo_adjoint_flow(ARGS[1], ARGS[2]; maxdim = length(ARGS)==3 ? parse(Int, ARGS[3]) : 16)
    end
end
