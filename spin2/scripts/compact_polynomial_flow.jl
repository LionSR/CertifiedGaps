import ThermalMPOGap
using JSON3, SHA
const MS = ThermalMPOGap.CertifiedMPO.Squaring
const MN = MS.MN

"Replace exact rational tokens, preserving every decoded value and other field."
function compact_tokens!(data)
    if data isa AbstractDict
        for key in keys(data)
            data[key] = compact_tokens!(data[key])
        end
    elseif data isa AbstractVector
        for index in eachindex(data)
            data[index] = compact_tokens!(data[index])
        end
    elseif data isa AbstractString && occursin(r"^[+-]?\d+//(?:[+-]?\d+|2\^\d+)$", data)
        value = MN.parse_fraction(data)
        token = MN.encode_fraction(value)
        @assert MN.parse_fraction(token) == value
        return token
    end
    return data
end

"Save a smaller copy of a polynomial flow; the standalone verifier must replay it."
function compact_polynomial_flow(source, output)
    @assert abspath(source) != abspath(output)
    @assert abspath(dirname(source)) == abspath(dirname(output))
    @assert !isfile(output)
    result = JSON3.read(read(source, String), Dict{String,Any})
    original_bytes, compact_bytes = 0, 0
    for row in result["runs"]
        input_path = joinpath(dirname(source), row["witness_file"])
        input_text = read(input_path, String)
        @assert bytes2hex(sha256(codeunits(input_text))) == row["witness_sha256"]
        data = JSON3.read(input_text, Dict{String,Any})
        compact_tokens!(data)
        output_path = splitext(output)[1]*"-step$(row["step"])-witness.json"
        @assert !isfile(output_path)
        open(output_path, "w") do io
            JSON3.write(io, data)
        end
        original_bytes += sizeof(input_text)
        compact_bytes += filesize(output_path)
        row["witness_file"] = basename(output_path)
        row["witness_sha256"] = bytes2hex(sha256(read(output_path)))
        println(
            "Compacted step $(row["step"]): $(sizeof(input_text)) -> $(filesize(output_path)) bytes",
        )
        flush(stdout)
    end
    result["reencoded_from"] = source
    result["reencoded_from_sha256"] = bytes2hex(sha256(read(source)))
    result["reencoding"] = Dict(
        "method"=>"lossless-dyadic-denominator-tokens",
        "original_witness_bytes"=>original_bytes,
        "compact_witness_bytes"=>compact_bytes,
        "independent_replay_required"=>true,
    )
    open(output, "w") do io
        JSON3.write(io, result)
    end
    return result
end

if abspath(PROGRAM_FILE) == @__FILE__
    @assert length(ARGS) == 2
    compact_polynomial_flow(ARGS...)
end
