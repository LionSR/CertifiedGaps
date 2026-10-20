include(joinpath(@__DIR__, "verify_polynomial_seed.jl"))
using SHA

"Condition an exact seed without introducing an initial approximation error."
function condition_polynomial_seed(source, output; bits = 128, precision = 384)
    started = time()
    @assert verify_seed(source) == 0
    data = JSON3.read(read(source, String), Dict{String,Any})
    @assert get(data, "exact_compression", false) && !get(data, "exact_conditioning", false)
    original = MN.decode_cores(data["candidate"])
    conditioned = PM.condition_right_exact(original; bits, precision)
    data["source_file"] = source
    data["source_sha256"] = bytes2hex(sha256(read(source)))
    data["unconditioned_candidate"] = data["candidate"]
    data["candidate"] = MN.encode_cores(conditioned.cores)
    data["exact_conditioning"] = true
    data["conditioning_bits"] = bits
    data["conditioning_transforms"] = PM.encode_transforms(conditioned.transforms)
    data["precision_bits"] = precision
    bounds = [
        MN.core_norm_upper(c, :right_hs; weights = PM.WEIGHTS, precision) for c in conditioned.cores
    ]
    data["right_hs_norm_bounds"] = string.(bounds)
    data["seconds"] = time()-started
    open(output, "w") do io
        JSON3.write(io, data)
    end
    @assert verify_seed(output) == 0
    println(
        "Exact conditioned seed: right HS bounds $(Float64.(bounds)); seconds=$(time()-started)",
    )
    return data
end

if abspath(PROGRAM_FILE) == @__FILE__
    @assert length(ARGS) == 2
    condition_polynomial_seed(ARGS...)
end
