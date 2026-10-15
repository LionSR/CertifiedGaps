include(joinpath(@__DIR__, "../polynomial/verify_polynomial_flow.jl"))
import CertifiedGaps
using CertifiedGaps.Results: write_result
const OE = CertifiedGaps.CertifiedMPO.OperatorError
const SB = CertifiedGaps.CertifiedMPO.SchattenBridge

"Replay an existing Gram flow and compute a backward Schatten trace bound."
function evaluate_schatten_flow(path; operator_files = Dict{Int,String}())
    verified = verify_flow(path)
    flow = JSON3.read(read(path, String))
    source = JSON3.read(read(String(flow.source_file), String))
    n = Int(source.length)
    @assert !isempty(flow.runs)
    @assert all(j -> 1 <= j <= length(flow.runs), keys(operator_files))
    kappas, hs_errors = MN.F[], MN.F[]
    operator_errors = Union{Nothing,MN.F}[]
    operator_sources = Any[]
    for (j, row) in enumerate(flow.runs)
        witness_path = joinpath(dirname(path), String(row.witness_file))
        proof = JSON3.read(read(witness_path, String))
        input, square = MS.decode_witness(proof)
        @assert square.adjoint_first
        kappa = MN.parse_fraction(proof.normalization)
        rho = MN.parse_fraction(proof.canonicalization_error_upper)
        push!(kappas, kappa)
        # Recanonicalization is controlled in HS norm; it also bounds every
        # Schatten p>=2 norm, including the operator norm.
        push!(hs_errors, square.error+kappa*rho)
        if haskey(operator_files, j)
            file = operator_files[j]
            text = read(file, String)
            encoded = JSON3.read(text)
            @assert encoded.source_sha256 == row.witness_sha256
            @assert bytes2hex(sha256(read(String(encoded.source_file)))) == row.witness_sha256
            operator = OE.decode(encoded)
            @assert OE.verify_operator_error(
                input,
                square,
                operator;
                residual_arithmetic = :interval,
            )
            push!(operator_errors, operator.error+kappa*rho)
            push!(
                operator_sources,
                Dict("step"=>j, "file"=>file, "sha256"=>bytes2hex(sha256(codeunits(text)))),
            )
        else
            push!(operator_errors, nothing)
        end
    end
    last_path = joinpath(dirname(path), String(last(flow.runs).witness_file))
    terminal = MN.parse_fraction(JSON3.read(read(last_path, String)).next_hs_norm_upper)
    initial_error = MN.parse_fraction(source.error_upper)
    chain = SB.gram_chain_upper(
        terminal,
        kappas,
        hs_errors;
        dimension = big(5)^n,
        operator_errors,
        initial_hs_error = initial_error,
    )
    @assert chain.initial_norm > 0
    RR = Nemo.ArbField(256)
    shift = MN.parse_fraction(flow.shift)
    # A_0 = Q/sqrt(5^N), Q >= exp(-delta*(H+N)),
    # p=2^(k+1), beta=delta*p. No approximate eigenvector enters this step.
    logtrace =
        -RR(shift-1)*RR(verified.beta)*n +
        RR(chain.power)*(log(RR(chain.initial_norm))+RR(n)*log(RR(5))/2)
    upper = PM.PG.rational_upper(logtrace)
    @assert RR(upper) > logtrace
    if n == 4
        spectrum = [
            (div(s*(s+1)-a*(a+1)-c*(c+1), 2), 2s+1) for a in 0:4 for c in 0:4 for
            s in abs(a-c):(a+c)
        ]
        exact =
            log(sum(deg*exp(-RR(verified.beta)*e) for (e, deg) in spectrum)) -
            RR(shift)*RR(verified.beta)*n
        @assert RR(upper) > exact
    end
    return (;
        upper,
        chain,
        operator_sources,
        beta = verified.beta,
        target_reached = verified.target_reached,
        length = n,
    )
end

function refine_schatten_flow(source, output; operator_files = Dict{Int,String}())
    @assert abspath(source) != abspath(output)
    checked = evaluate_schatten_flow(source; operator_files)
    source_text = read(source, String)
    original = JSON3.read(source_text)
    baseline = MN.parse_fraction(
        get(
            original,
            :backward_shifted_log_trace_upper,
            last(original.runs).shifted_log_trace_upper,
        ),
    )
    result = Dict(
        "status"=>"backward-Schatten-Gram-trace-certificate",
        "source_file"=>source,
        "source_sha256"=>bytes2hex(sha256(codeunits(source_text))),
        "length"=>checked.length,
        "steps"=>length(original.runs),
        "beta"=>string(checked.beta),
        "moment_power"=>string(checked.chain.power),
        "candidate_norm_bounds"=>string.(checked.chain.candidate_norms),
        "dimension_root_bounds"=>string.(checked.chain.dimension_roots),
        "local_schatten_error_bounds"=>string.(checked.chain.local_errors),
        "initial_schatten_norm_upper"=>string(checked.chain.initial_norm),
        "operator_sources"=>checked.operator_sources,
        "shifted_log_trace_upper"=>string(checked.upper),
        "shifted_log_trace_upper_display"=>Float64(checked.upper),
        "previous_shifted_log_trace_upper"=>string(baseline),
        "log_bound_change_display"=>Float64(checked.upper-baseline),
        "current_thermal_trace_upper_certified"=>true,
        "target_beta_reached"=>checked.target_reached,
        "uniform_gap_certified"=>false,
    )
    write_result(output, result)
    println(
        "Schatten trace bound at N=$(checked.length), beta=$(Float64(checked.beta)): log upper $(Float64(checked.upper)); change $(Float64(checked.upper-baseline))",
    )
    return result
end

"Independently replay all source inputs and the saved Schatten conclusion."
function verify_schatten_flow(path)
    data = JSON3.read(read(path, String))
    source = String(data.source_file)
    @assert bytes2hex(sha256(read(source))) == data.source_sha256
    operator_files = Dict{Int,String}()
    for row in data.operator_sources
        @assert !haskey(operator_files, Int(row.step))
        @assert bytes2hex(sha256(read(String(row.file)))) == row.sha256
        operator_files[Int(row.step)] = String(row.file)
    end
    checked = evaluate_schatten_flow(source; operator_files)
    @assert Int(data.length) == checked.length
    @assert Int(data.steps) >= 0 && big(2)^(Int(data.steps)+1) == checked.chain.power
    @assert MN.parse_fraction(data.beta) == checked.beta
    @assert parse(BigInt, data.moment_power) == checked.chain.power
    @assert MN.parse_fraction.(data.candidate_norm_bounds) == checked.chain.candidate_norms
    @assert MN.parse_fraction.(data.dimension_root_bounds) == checked.chain.dimension_roots
    @assert MN.parse_fraction.(data.local_schatten_error_bounds) == checked.chain.local_errors
    @assert MN.parse_fraction(data.initial_schatten_norm_upper) == checked.chain.initial_norm
    @assert MN.parse_fraction(data.shifted_log_trace_upper) >= checked.upper
    @assert data.target_beta_reached == checked.target_reached && !data.uniform_gap_certified
    println("Verified backward Schatten thermal trace certificate: $path")
    return checked
end

if abspath(PROGRAM_FILE) == @__FILE__
    @assert length(ARGS) == 2
    refine_schatten_flow(ARGS...)
end
