# Certify squaring/compression error for quantized thermal candidates.
# Their initial distance from the physical Gibbs amplitude is not certified here.
include(joinpath(@__DIR__, "mpo_gauge_probe.jl"))
import ThermalMPOGap
using ThermalMPOGap.Results: write_result
const MS = ThermalMPOGap.CertifiedMPO.Squaring

function squaring_main(
    output = joinpath(@__DIR__, "../results/mpo-squaring-probe.json");
    lengths = (4, 8),
    maxdim = 8,
    beta = 0.125,
    optimize_norm_gauge = false,
    blocksize = 1,
)
    result = Dict{String,Any}(
        "status"=>"certified-squaring-of-exact-candidate-MPO",
        "complete"=>false,
        "thermal_trace_certified"=>false,
        "uniform_gap_certified"=>false,
        "runs"=>Any[],
        "runner_sha256"=>bytes2hex(sha256(read(@__FILE__))),
        "squaring_module_sha256"=>bytes2hex(
            sha256(read(joinpath(@__DIR__, "../../julia/src/CertifiedMPO/Squaring.jl"))),
        ),
    )
    for n in lengths
        sites = siteinds("S=2", n; conserve_qns = true)
        amplitude = MPO(sites, "Id")
        amplitude /= norm(amplitude)
        half = [exp(-beta/4*Spin2Gap.bond_operator(sites, j, k)) for (j, k) in bonds(n)]
        amplitude = apply(vcat(half, reverse(half)), amplitude; maxdim, cutoff = 1e-12)
        amplitude /= norm(amplitude)
        cores = MN.quantize(extract_cores(orthogonalize!(amplitude, 1), sites))
        started = time()
        gauges =
            optimize_norm_gauge ?
            MN.rational_gauges(cores, optimize_gauge([Float64.(c) for c in cores]).x) : nothing
        certificate = MS.square_candidate(cores; maxdim, gauges, blocksize)
        @assert MS.verify_square_certificate(cores, certificate)
        norm2 = MS.hs_squared(certificate.cores)
        norm_lower = norm2/MN.relative_upper_sqrt(norm2)
        relative_error = certificate.error/norm_lower
        row = Dict{String,Any}(
            "length"=>n,
            "initial_beta"=>beta,
            "optimized_first_copy_gauge"=>optimize_norm_gauge,
            "norm_blocksize"=>blocksize,
            "maximum_bond_dimension"=>maxdim,
            "compression_error_hs_upper"=>string(certificate.error),
            "relative_to_candidate_hs_upper"=>Float64(relative_error),
            "squaring_error_certified"=>true,
            "independent_replay_passed"=>true,
            "local_error_terms"=>[Float64(r.error_term) for r in certificate.rows],
            "seconds"=>time()-started,
        )
        if n == 4
            B, C = dense_operator(cores), dense_operator(certificate.cores)
            row["independent_dense_error"] = norm(B*B-C)
            row["independent_dense_relative_error"] = norm(B*B-C)/norm(C)
        end
        proof = Dict(
            "input"=>MN.encode_cores(cores),
            "output"=>MN.encode_cores(certificate.cores),
            "bridges"=>[
                Dict("dimensions"=>collect(size(B)), "entries"=>string.(vec(B))) for
                B in certificate.bridges
            ],
            "residual_bounds"=>[string(r.residual_bound) for r in certificate.rows],
            "prefix_bounds"=>[string(r.local_prefix_bound) for r in certificate.rows],
            "operator_bounds"=>string.(certificate.tail.local_operator_bounds),
            "hs_bounds"=>string.(certificate.tail.local_hs_bounds),
            "tail_frobenius_capped"=>certificate.tail.frobenius_capped,
            "block_bounds"=>string.(certificate.tail.block_bounds),
            "error"=>string(certificate.error),
            "precision_bits"=>192,
            "thermal_trace_certified"=>false,
            "uniform_gap_certified"=>false,
        )
        if !isnothing(gauges)
            proof["gauges"] = [string.(g) for g in gauges]
        end
        proofpath = splitext(output)[1]*"-n$(n)-witness.json"
        open(proofpath, "w") do io
            JSON3.write(io, proof)
        end
        row["witness_file"], row["witness_sha256"] =
            basename(proofpath), bytes2hex(sha256(read(proofpath)))
        push!(result["runs"], row)
        write_result(output, result)
        println(
            "squaring n=$n D=$maxdim relative error <= $(Float64(relative_error)) seconds=$(time()-started)",
        )
        flush(stdout)
    end
    result["complete"] = true
    write_result(output, result)
end

if abspath(PROGRAM_FILE) == @__FILE__
    squaring_main(
        isempty(ARGS) ? joinpath(@__DIR__, "../results/mpo-squaring-probe.json") : only(ARGS),
    )
end
