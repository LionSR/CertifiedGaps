import ThermalMPOGap
using JSON3
const MS = ThermalMPOGap.CertifiedMPO.Squaring
const MN = MS.MN

function verify_saved_squaring(path)
    data = JSON3.read(read(path, String))
    input, certificate = MS.decode_witness(data)
    @assert MS.verify_square_certificate(input, certificate; precision = Int(data.precision_bits))
    println("Verified MPO squaring error <= $(Float64(certificate.error)): $path")
    return certificate.error
end

if abspath(PROGRAM_FILE) == @__FILE__
    for path in ARGS
        verify_saved_squaring(path)
    end
end
