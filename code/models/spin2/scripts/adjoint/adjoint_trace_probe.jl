# Dense four-site reference for the full multi-step adjoint correction.
include(joinpath(@__DIR__, "../compression/repeated_l2_probe.jl"))
import CertifiedGaps
using CertifiedGaps.Results: write_result
const AT = CertifiedGaps.CertifiedMPO.AdjointTrace

grid_upper(x::F; bits = 80) = cld(numerator(x)*big(2)^bits, denominator(x))//big(2)^bits

function jordan_interval(B, G)
    RR=Nemo.ArbField(192)
    lower, upper=zero(B), zero(B)
    for states in SECTORS
        b=Nemo.matrix(RR, [B[location(r, c)...] for r in states, c in states])
        g=Nemo.matrix(RR, [G[location(r, c)...] for r in states, c in states])
        product=b*g+g*b
        for (i, r) in enumerate(states), (j, c) in enumerate(states)
            lower[location(r, c)...], upper[location(r, c)...]=PB.rational_bracket(product[i, j])
        end
    end
    return hermitize(lower), hermitize(upper)
end

function trace_pair_ball(A, B)
    RR=Nemo.ArbField(192)
    total=RR(0)
    for states in SECTORS
        a=Nemo.matrix(RR, [A[location(r, c)...] for r in states, c in states])
        b=Nemo.matrix(RR, [B[location(r, c)...] for r in states, c in states])
        total+=Nemo.tr(a*b)
    end
    return total
end

function signed_center(initial, B, G)
    RR=Nemo.ArbField(192)
    center=RR(norm2(B[end]))+2trace_pair_ball(G[1], initial-B[1])
    for j in 1:(length(B)-1)
        for states in SECTORS
            b=Nemo.matrix(RR, [B[j][location(r, c)...] for r in states, c in states])
            next=Nemo.matrix(RR, [B[j+1][location(r, c)...] for r in states, c in states])
            g=Nemo.matrix(RR, [G[j+1][location(r, c)...] for r in states, c in states])
            center+=2Nemo.tr(g*(b*b-next))
        end
    end
    return center
end

encode_matrix(M) = Dict("entries"=>[[i, string(M[i])] for i in eachindex(M) if !iszero(M[i])])
fraction(s) = begin
    a, b=split(String(s), "//")
    parse(BigInt, a)//parse(BigInt, b)
end
function decode_matrix(data)
    M=zeros(F, 625, 625)
    seen=Set{Int}()
    for row in data["entries"]
        i=Int(row[1])
        @assert 1<=i<=length(M) && !(i in seen)
        push!(seen, i)
        M[i]=fraction(row[2])
    end
    return M
end
encode_norm(p) = Dict("upper"=>string(p.upper), "vectors"=>[string.(v) for v in p.vectors])
decode_norm(p) = (; upper = fraction(p["upper"]), vectors = [fraction.(v) for v in p["vectors"]])

function evaluate_flow(initial, initial_rounding, B, G, errors, tau, dual_norms)
    center=signed_center(initial, B, G)
    lo, hi=PB.rational_bracket(center)
    initial_uncertainty=2initial_rounding*sqrt_upper(norm2(G[1]))
    return AT.trace_interval(
        lo,
        hi;
        errors,
        dual_defects = tau,
        dual_operator_bounds = [p.upper for p in dual_norms],
        normalizations = ones(F, length(B)-1),
        initial_uncertainty,
    )
end

function run_adjoint(
    output;
    rankcaps = (32, 64),
    beta_target = 4,
    doublings = 5,
    dual_rankcap = nothing,
)
    result=Dict{String,Any}(
        "status"=>"four-site-multistep-adjoint-trace-certificate",
        "complete"=>false,
        "uniform_gap_certified"=>false,
        "large_chain_trace_certified"=>false,
        "runs"=>Any[],
        "beta_target"=>beta_target,
        "doublings"=>doublings,
        "dual_rankcap"=>dual_rankcap,
    )
    lower, upper, _=amplitude_interval(beta_target//2^doublings, 2)
    initial=hermitize((lower+upper)/2)
    initial_rounding=sqrt_upper(norm2((upper-lower)/2))
    reference=last(amplitude_interval(beta_target, 2^(doublings+1)))
    for rankcap in rankcaps
        started=time()
        println("Starting adjoint trace cap=$rankcap, beta=$beta_target, steps=$doublings")
        flush(stdout)
        initial_candidate, rank=svd_candidate(initial, rankcap)
        B=[initial_candidate]
        errors=F[grid_upper(initial_rounding+sqrt_upper(norm2(initial-initial_candidate)))]
        primal_norms=Any[]
        for j in 1:doublings
            norm_witness=operator_norm_upper(B[end]; details = true)
            @assert verify_operator_norm_witness(B[end], norm_witness)
            push!(primal_norms, norm_witness)
            lower, upper=squared_interval(B[end])
            midpoint=hermitize((lower+upper)/2)
            next, rank=svd_candidate(midpoint, rankcap)
            delta=sqrt_upper(norm2(midpoint-next))+sqrt_upper(norm2((upper-lower)/2))
            push!(errors, grid_upper(delta+errors[end]*(2norm_witness.upper+errors[end])))
            push!(B, next)
        end
        println("Forward flow complete: error=$(Float64(errors[end]))")
        flush(stdout)
        G=Vector{Matrix{F}}(undef, doublings+1)
        G[end]=B[end]
        tau=zeros(F, doublings)
        dual_norms=Vector{Any}(undef, doublings)
        for j in doublings:-1:1
            dual_norms[j]=operator_norm_upper(G[j+1]; details = true)
            @assert verify_operator_norm_witness(G[j+1], dual_norms[j])
            lower, upper=jordan_interval(B[j], G[j+1])
            midpoint=hermitize((lower+upper)/2)
            G[j]=isnothing(dual_rankcap) ? midpoint : first(svd_candidate(midpoint, dual_rankcap))
            tau[j]=grid_upper(sqrt_upper(norm2((upper-lower)/2))+sqrt_upper(norm2(midpoint-G[j])))
        end
        bounds=evaluate_flow(initial, initial_rounding, B, G, errors, tau, dual_norms)
        _, norm_hi=PB.rational_bracket(sqrt(Nemo.ArbField(192)(norm2(B[end]))))
        triangle=(norm_hi+last(errors))^2
        @assert reference < bounds.upper && reference > bounds.lower && reference < triangle
        proof=Dict(
            "beta_target"=>beta_target,
            "doublings"=>doublings,
            "primal"=>encode_matrix.(B),
            "dual"=>encode_matrix.(G),
            "primal_errors"=>string.(errors),
            "dual_defects"=>string.(tau),
            "primal_norms"=>encode_norm.(primal_norms),
            "dual_norms"=>encode_norm.(dual_norms),
            "trace_lower"=>string(bounds.lower),
            "trace_upper"=>string(bounds.upper),
            "uniform_gap_certified"=>false,
        )
        proofpath=splitext(output)[1]*"-cap$(rankcap)-witness.json"
        open(io->JSON3.write(io, proof), proofpath, "w")
        row=Dict(
            "rankcap_per_charge"=>rankcap,
            "final_rank_upper"=>rank,
            "triangle_over_reference_display"=>Float64(triangle)/Float64(reference),
            "corrected_over_reference_display"=>Float64(bounds.upper)/Float64(reference),
            "corrected_lower_over_reference_display"=>Float64(bounds.lower)/Float64(reference),
            "quadratic_remainder_display"=>Float64(bounds.quadratic),
            "primal_dual_remainder_display"=>Float64(bounds.cross),
            "terminal_remainder_display"=>Float64(bounds.terminal),
            "reference_trace_ball"=>string(reference),
            "witness_file"=>basename(proofpath),
            "witness_sha256"=>bytes2hex(sha256(read(proofpath))),
            "seconds"=>time()-started,
        )
        push!(result["runs"], row)
        write_result(output, result)
        println(
            "cap=$rankcap triangle/reference=$(row["triangle_over_reference_display"]) corrected/reference=$(row["corrected_over_reference_display"]) seconds=$(time()-started)",
        )
        flush(stdout)
    end
    result["complete"]=true
    write_result(output, result)
end

function verify_adjoint(path)
    data=JSON3.read(read(path, String))
    steps, beta=Int(data.doublings), Int(data.beta_target)
    B, G=decode_matrix.(data.primal), decode_matrix.(data.dual)
    errors, tau=fraction.(data.primal_errors), fraction.(data.dual_defects)
    primal_norms, dual_norms=decode_norm.(data.primal_norms), decode_norm.(data.dual_norms)
    @assert length(B)==length(G)==length(errors)==steps+1 && length(tau)==steps
    @assert length(primal_norms)==length(dual_norms)==steps
    @assert G[end]==B[end] && all(hermitize(x)==x for x in vcat(B, G))
    lower, upper, _=amplitude_interval(beta//2^steps, 2)
    initial=hermitize((lower+upper)/2)
    rounding=sqrt_upper(norm2((upper-lower)/2))
    @assert errors[1]>=rounding+sqrt_upper(norm2(initial-B[1]))
    for j in 1:steps
        @assert verify_operator_norm_witness(B[j], primal_norms[j])
        @assert verify_operator_norm_witness(G[j+1], dual_norms[j])
        lower, upper=squared_interval(B[j])
        midpoint=hermitize((lower+upper)/2)
        delta=sqrt_upper(norm2(midpoint-B[j+1]))+sqrt_upper(norm2((upper-lower)/2))
        @assert errors[j+1]>=delta+errors[j]*(2primal_norms[j].upper+errors[j])
        lower, upper=jordan_interval(B[j], G[j+1])
        midpoint=hermitize((lower+upper)/2)
        @assert tau[j]>=sqrt_upper(norm2(midpoint-G[j]))+sqrt_upper(norm2((upper-lower)/2))
    end
    bounds=evaluate_flow(initial, rounding, B, G, errors, tau, dual_norms)
    @assert fraction(data.trace_upper)>=bounds.upper && fraction(data.trace_lower)<=bounds.lower
    reference=last(amplitude_interval(beta, 2^(steps+1)))
    @assert reference < fraction(data.trace_upper) && reference > fraction(data.trace_lower)
    println("Verified all $steps squarings and adjoint corrections: $path")
    return true
end

if abspath(PROGRAM_FILE)==@__FILE__
    if first(ARGS)=="--verify"
        foreach(verify_adjoint, ARGS[2:end])
    else
        run_adjoint(only(ARGS))
    end
end
