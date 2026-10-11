# Reuse saved exact inputs and gauges; increase the output rank without rerunning ITensors.
include(joinpath(@__DIR__, "verify_mpo_squaring.jl"))
using SHA
using CertifiedGaps.Results: write_result

function refine_squaring(source, output; rankcaps = (24, 32), bits = 48)
    data = JSON3.read(read(source, String))
    input = MN.decode_cores(data.input)
    gauges = haskey(data, :gauges) ? MN.decode_gauges(data.gauges) : nothing
    result = Dict{String,Any}(
        "status"=>"certified-squaring-rank-refinement",
        "complete"=>false,
        "source_sha256"=>bytes2hex(sha256(read(source))),
        "source_file"=>source,
        "thermal_trace_certified"=>false,
        "uniform_gap_certified"=>false,
        "runner_sha256"=>bytes2hex(sha256(read(@__FILE__))),
        "runs"=>Any[],
    )
    for maxdim in rankcaps
        started = time()
        certificate = MS.square_candidate(input; maxdim, bits, gauges, blocksize = 2)
        @assert MS.verify_square_certificate(input, certificate)
        norm2 = MS.hs_squared(certificate.cores)
        relative = certificate.error*MN.relative_upper_sqrt(norm2)/norm2
        proof = MS.witness_encoding(input, certificate)
        proofpath = splitext(output)[1]*"-D$(maxdim)-witness.json"
        open(proofpath, "w") do io
            JSON3.write(io, proof)
        end
        row = Dict(
            "length"=>length(input),
            "output_bond_dimension"=>maxdim,
            "rounding_bits"=>bits,
            "error_hs_upper"=>string(certificate.error),
            "relative_to_candidate_hs_upper"=>Float64(relative),
            "squaring_error_certified"=>true,
            "independent_replay_passed"=>true,
            "witness_file"=>basename(proofpath),
            "witness_sha256"=>bytes2hex(sha256(read(proofpath))),
            "seconds"=>time()-started,
        )
        push!(result["runs"], row)
        write_result(output, result)
        println(
            "refinement n=$(length(input)) output D=$maxdim relative error <= $(Float64(relative)) seconds=$(time()-started)",
        )
        flush(stdout)
    end
    result["complete"] = true
    write_result(output, result)
end

if abspath(PROGRAM_FILE) == @__FILE__
    @assert length(ARGS) == 2
    refine_squaring(ARGS...)
end
