import CertifiedGaps
using JSON3, SHA
const OE = CertifiedGaps.CertifiedMPO.OperatorError

function verify_operator_witness(path; residual_arithmetic = :interval)
    data = JSON3.read(read(path, String))
    source = String(data.source_file)
    @assert bytes2hex(sha256(read(source))) == data.source_sha256
    input, square = OE.MS.decode_witness(JSON3.read(read(source, String)))
    proof = OE.decode(data)
    progress = row -> begin
        if row.site % 16 == 0 || row.site == row.total-1
            println("Replay operator residual: site $(row.site)/$(row.total-1)")
            flush(stdout)
        end
    end
    @assert OE.verify_operator_error(
        input,
        square,
        proof;
        residual_arithmetic,
        precision = get(data, :precision_bits, 192),
        progress,
    )
    println("Verified operator compression error <= $(Float64(proof.error)): $path")
end

if abspath(PROGRAM_FILE) == @__FILE__
    foreach(verify_operator_witness, ARGS)
end
