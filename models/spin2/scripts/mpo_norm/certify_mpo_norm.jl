# Certify the norm of saved exact rational MPO candidates; not a thermal error certificate.
include(joinpath(@__DIR__, "mpo_gauge_probe.jl"))
using CertifiedGaps.Results: write_result

function certify_main(
    output = joinpath(@__DIR__, "../../results/mpo-local-norm-certificate.json");
    lengths = (4, 8),
    target_steps = (1, 8),
    iterations = 200,
    conserve_qns = true,
)
    result = Dict{String,Any}(
        "status"=>"exact-candidate-MPO-norm-certificate",
        "complete"=>false,
        "thermal_trace_certified"=>false,
        "uniform_gap_certified"=>false,
        "conserve_qns"=>conserve_qns,
        "runs"=>Any[],
        "runner_sha256"=>bytes2hex(sha256(read(@__FILE__))),
        "norm_module_sha256"=>bytes2hex(
            sha256(read(joinpath(@__DIR__, "../../../../julia/src/CertifiedMPO/Norms.jl"))),
        ),
    )
    for n in lengths
        sites = siteinds("S=2", n; conserve_qns)
        amplitude = MPO(sites, "Id")
        amplitude /= norm(amplitude)
        half = [exp(-0.125/4*Spin2Gap.bond_operator(sites, j, k)) for (j, k) in bonds(n)]
        gates = vcat(half, reverse(half))
        for step in 1:maximum(target_steps)
            amplitude = apply(gates, amplitude; maxdim = 64, cutoff = 1e-12)
            amplitude /= norm(amplitude)
            step in target_steps || continue
            started = time()
            numerical = extract_cores(orthogonalize!(deepcopy(amplitude), div(n, 2)), sites)
            cores = MN.quantize(numerical)
            best = optimize_gauge(numerical; iterations)
            gauges = MN.rational_gauges(cores, best.x)
            certificate = MN.local_product_witness(cores, gauges)
            row = Dict{String,Any}(
                "length"=>n,
                "beta"=>step/8,
                "maximum_bond_dimension"=>maxlinkdim(amplitude),
                "local_norm_product_initial"=>exp(best.initial),
                "local_norm_product_numerical"=>exp(best.final),
                "candidate_operator_norm_upper"=>string(certificate.upper),
                "candidate_operator_norm_upper_display"=>Float64(certificate.upper),
                "gradient_norm"=>best.gradient_norm,
                "norm_bound_certified"=>true,
                "precision_bits"=>192,
            )
            if n == 4
                row["independent_dense_candidate_norm"] = opnorm(dense_operator(cores))
            end
            witness = Dict(
                "status"=>"exact-rational-candidate-MPO-norm-witness",
                "precision_bits"=>192,
                "cores"=>MN.encode_cores(cores),
                "gauges"=>[string.(g) for g in gauges],
                "local_bounds"=>string.(certificate.bounds),
                "upper"=>string(certificate.upper),
                "thermal_trace_certified"=>false,
                "uniform_gap_certified"=>false,
            )
            witness_path = splitext(output)[1]*"-n$(n)-step$(step)-witness.json"
            open(witness_path, "w") do io
                JSON3.write(io, witness)
            end
            row["witness_file"] = basename(witness_path)
            row["witness_sha256"] = bytes2hex(sha256(read(witness_path)))
            row["seconds"] = time()-started
            push!(result["runs"], row)
            write_result(output, result)
            println(
                "local norm n=$n beta=$(step/8) initial=$(exp(best.initial)) certified=$(Float64(certificate.upper)) seconds=$(time()-started)",
            )
            flush(stdout)
        end
    end
    result["complete"] = true
    result["completed_utc"] = string(now(UTC))
    write_result(output, result)
end

if abspath(PROGRAM_FILE) == @__FILE__
    certify_main(
        isempty(ARGS) ? joinpath(@__DIR__, "../../results/mpo-local-norm-certificate.json") :
        only(ARGS),
    )
end
