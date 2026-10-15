using Spin2Gap, ITensors, ITensorMPS, LinearAlgebra, JSON3, SHA
const PM = Spin2Gap.PolynomialMPO
const MS = PM.MS
const MN = PM.MN
BLAS.set_num_threads(1)

function candidate_mpo(cores; maxdim = 16, bits = 48)
    n = length(cores)
    sites = [Index(5, "Site,n=$j") for j in 1:n]
    links = [Index(size(cores[j], 4), "Link,l=$j") for j in 1:(n-1)]
    physical = sqrt.([Float64(PM.WEIGHTS[t]/PM.WEIGHTS[s]) for s in 1:5, t in 1:5])
    tensors = ITensor[]
    for (j, core) in enumerate(cores)
        indices = vcat(
            j==1 ? Index[] : [links[j-1]],
            [prime(sites[j]), dag(sites[j])],
            j==n ? Index[] : [links[j]],
        )
        A = Float64.(core) .* reshape(physical, 1, 5, 5, 1)
        push!(tensors, ITensor(reshape(A, dim.(indices)...), indices...))
    end
    A = MPO(tensors)
    truncate!(A; maxdim, cutoff = 0.0)
    orthogonalize!(A, 1)
    result = Array{MN.F,4}[]
    for j in 1:n
        left = j==1 ? Index[] : [linkind(A, j-1)]
        right = j==n ? Index[] : [linkind(A, j)]
        raw = Array(A[j], left..., prime(sites[j]), dag(sites[j]), right...)
        core = reshape(
            raw,
            isempty(left) ? 1 : dim(only(left)),
            5,
            5,
            isempty(right) ? 1 : dim(only(right)),
        )
        push!(result, MS.relative_quantize(core ./ reshape(physical, 1, 5, 5, 1); bits))
    end
    return result
end

function seed_probe(
    output = joinpath(@__DIR__, "../../results/polynomial-seed-probe.json");
    lengths = (4, 16),
    maxdim = 16,
    power = 28,
    beta = 64,
)
    result = Dict{String,Any}(
        "status"=>"initial-polynomial-MPO-error-certificate",
        "complete"=>false,
        "uniform_gap_certified"=>false,
        "thermal_trace_certified"=>false,
        "runs"=>Any[],
        "runner_sha256"=>bytes2hex(sha256(read(@__FILE__))),
    )
    for n in lengths
        started = time()
        seed = PM.PG.seed(n, beta, power)
        exact = PM.normalized_polynomial(seed)
        candidate = candidate_mpo(exact; maxdim)
        check = PM.initial_error(exact, candidate)
        norm_lower = PM.PG.rational_upper(-sqrt(check.cc))*(-1)
        @assert norm_lower > 0
        proof = Dict(
            "length"=>n,
            "beta"=>string(seed.beta),
            "power"=>power,
            "remainder"=>string(seed.remainder),
            "candidate"=>MN.encode_cores(candidate),
            "error_upper"=>string(check.error),
            "precision_bits"=>192,
            "thermal_trace_certified"=>false,
            "uniform_gap_certified"=>false,
        )
        proofpath = splitext(output)[1]*"-n$(n)-witness.json"
        open(proofpath, "w") do io
            JSON3.write(io, proof)
        end
        row = Dict(
            "length"=>n,
            "maxdim"=>maxdim,
            "power"=>power,
            "beta"=>beta,
            "initial_error_certified"=>true,
            "initial_error_upper"=>string(check.error),
            "relative_initial_error_upper"=>Float64(check.error/norm_lower),
            "candidate_norm_interval"=>string(sqrt(check.cc)),
            "witness_file"=>basename(proofpath),
            "witness_sha256"=>bytes2hex(sha256(read(proofpath))),
            "seconds"=>time()-started,
        )
        push!(result["runs"], row)
        write_result(output, result)
        println(
            "polynomial seed n=$n D=$maxdim initial relative error <= $(Float64(check.error/norm_lower)) seconds=$(time()-started)",
        )
        flush(stdout)
    end
    result["complete"] = true
    write_result(output, result)
end

if abspath(PROGRAM_FILE) == @__FILE__
    seed_probe(
        isempty(ARGS) ? joinpath(@__DIR__, "../../results/polynomial-seed-probe.json") : only(ARGS),
    )
end
