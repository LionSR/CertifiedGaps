import CertifiedGaps
using Spin2Gap, ITensors, ITensorMPS, LinearAlgebra, JSON3, SHA, Dates
const MN = CertifiedGaps.CertifiedMPO.Norms
BLAS.set_num_threads(1)

function extract_cores(A, sites)
    result = Array{Float64,4}[]
    for j in eachindex(A)
        left = j == 1 ? Index[] : [linkind(A, j-1)]
        right = j == length(A) ? Index[] : [linkind(A, j)]
        raw = Array(A[j], left..., prime(sites[j]), dag(sites[j]), right...)
        @assert all(isreal, raw)
        push!(
            result,
            reshape(
                real.(raw),
                isempty(left) ? 1 : dim(only(left)),
                5,
                5,
                isempty(right) ? 1 : dim(only(right)),
            ),
        )
    end
    return result
end

function dense_operator(cores)
    partial = [ones(1, 1)]
    for core in cores
        l, d, _, r = size(core)
        next = [zeros(d*size(first(partial), 1), d*size(first(partial), 2)) for _ in 1:r]
        for a in 1:r, b in 1:l
            next[a] += kron(partial[b], Float64.(core[b, :, :, a]))
        end
        partial = next
    end
    return only(partial)
end

function main(
    output = joinpath(@__DIR__, "../results/mpo-norm-probe.json");
    conserve_qns = true,
    lengths = (4, 8),
    target_steps = (1, 8),
)
    result = Dict{String,Any}(
        "status"=>"norm-witness-for-quantized-MPO-candidates-only",
        "complete"=>false,
        "uniform_gap_certified"=>false,
        "thermal_trace_certified"=>false,
        "conserve_qns"=>conserve_qns,
        "lengths"=>collect(lengths),
        "target_steps"=>collect(target_steps),
        "runs"=>Any[],
        "runner_sha256"=>bytes2hex(sha256(read(@__FILE__))),
        "norm_module_sha256"=>bytes2hex(
            sha256(read(joinpath(@__DIR__, "../../../julia/src/CertifiedMPO/Norms.jl"))),
        ),
        "ITensors"=>string(pkgversion(ITensors)),
        "ITensorMPS"=>string(pkgversion(ITensorMPS)),
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
            for center in unique((1, div(n, 2), n))
                started = time()
                candidate = orthogonalize!(deepcopy(amplitude), center)
                cores = MN.quantize(extract_cores(candidate, sites))
                certificate = MN.witness(cores)
                row = Dict{String,Any}(
                    "length"=>n,
                    "beta"=>step/8,
                    "gauge_center"=>center,
                    "maximum_bond_dimension"=>maxlinkdim(candidate),
                    "candidate_operator_norm_upper"=>string(certificate.upper),
                    "candidate_operator_norm_upper_display"=>Float64(certificate.upper),
                    "quantization_bits"=>48,
                    "norm_rounding_bits"=>80,
                    "numerical_candidate_hs_norm"=>norm(candidate),
                    "seconds"=>time()-started,
                )
                if n == 4
                    physical = dense_operator(cores)
                    reference = opnorm(physical)
                    @assert Float64(certificate.upper) >= reference-1e-12
                    row["independent_dense_candidate_norm"] = reference
                    row["upper_over_numerical_reference"] = Float64(certificate.upper)/reference
                end
                push!(result["runs"], row)
                write_result(output, result)
                println(
                    "MPO norm n=$n beta=$(step/8) center=$center upper=$(Float64(certificate.upper)) seconds=$(time()-started)",
                )
                flush(stdout)
            end
        end
    end
    result["complete"] = true
    result["completed_utc"] = string(now(UTC))
    write_result(output, result)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(isempty(ARGS) ? joinpath(@__DIR__, "../results/mpo-norm-probe.json") : only(ARGS))
end
