module Spin1Certificates

using Dates, JSON3, LinearAlgebra, SHA
using CertifiedGaps.ExactArithmetic: RationalBig, frac, parse_rational
using CertifiedGaps.Results: parse_exact_json, read_exact_json, source_hashes, write_result
using CertifiedGaps.Recurrence: squaring_defect, bootstrap_coefficient

const SOURCE =
    normpath(joinpath(@__DIR__, "../../../references/periodic-gap-source/verification/computations"))
const MANIFEST_SHA = "2d2461e8f87ef9c3af62ce226412b106fcfb34f01d5139683a90708859e2dc92"
const UPSTREAM_COMMIT = "adc7f1241b42e322a6451854ab7e4b4c146bf78a"
# The core keeps its exact integer parser private; this is the same three-method refusal.
integer(x::AbstractString) = parse(BigInt, x)
integer(x::Integer) = BigInt(x)
integer(x) = error("Refusing inexact integer input: $(typeof(x))")

function verify_source()
    path = joinpath(SOURCE, "source-manifest.json")
    @assert bytes2hex(sha256(read(path))) == MANIFEST_SHA
    manifest = read_exact_json(path)
    count = 0
    for group in ("files", "supplementary_certificate_files")
        for (name, expected) in manifest[group]
            @assert bytes2hex(sha256(read(joinpath(SOURCE, name)))) == expected "Source changed: $name"
            count += 1
        end
    end
    @assert bytes2hex(sha256(read(joinpath(SOURCE, "../verify.py")))) ==
            manifest["verify_py_sha256"]
    return Dict("manifest_sha256" => MANIFEST_SHA, "verified_files" => count + 1)
end

include("thermal.jl")
include("thermal_kernels.jl")
include("arithmetic.jl")
include("trials.jl")
include("trials_su2.jl")
include("scalars.jl")
include("audit.jl")
include("itensors_checks.jl")

end
