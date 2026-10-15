isdefined(@__MODULE__, :TransferCertificate) || include("src/TransferCertificate.jl")
using CertifiedGaps.Results
include("bridge.jl")

"Replay all new certificates, checking input, implementation and shared-core fingerprints."
function check_transfer(path = joinpath(@__DIR__, "certificate.json"))
    tc = TransferCertificate
    data = parse_exact_json(read(path, String))
    data["format"] == "spin1_largeD transfer and Weyl certificate v1" || error("format")
    Dict(data["source_hashes"]) == tc.proof_hashes() || error("proof source hash mismatch")
    Dict(data["input_hashes"]) == tc.input_hashes() || error("spectral input hash mismatch")
    keep = parse(Int, data["partition_keep"])
    keep == 8192 || error("unsupported partition cutoff")
    spectra = Dict(D => (tc.reference_spectrum(D, 6), tc.reference_spectrum(D, 12))
        for D in (4//1, 1603//500))
    isempty(data["cases"]) && error("empty certificate")
    for record in data["cases"]
        a, b = spectra[parse(tc.Q, record["reference_D"])]
        row = tc.check_record(record, a, b; keep)
        println(record["id"], ": D in ", record["D_interval"], ", even L >= ",
            row.n0, ", gap >= ", Float64(parse(tc.Q, record["bound"])))
    end
    lifted = only(filter(c -> c["id"] == "space-lift-1-D4//1", data["cases"]))
    bridge = parse_exact_json(read(joinpath(@__DIR__, "results", "finite_bridge.json"), String))
    check_bridge(bridge, parse(tc.Q, lifted["bound"])) || error("finite bridge failed")
    println("D=4: the stronger n0=12 bound also holds for every even L >= 6 sites")
    return data
end

if abspath(PROGRAM_FILE) == @__FILE__
    check_transfer()
end
