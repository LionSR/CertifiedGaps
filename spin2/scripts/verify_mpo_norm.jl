# Replay saved bounds without ITensors, a floating eigensolver, or gauge optimization.
import ThermalMPOGap
using JSON3, SHA
const MN = ThermalMPOGap.CertifiedMPO.Norms

function verify_saved_norm(path)
    data = JSON3.read(read(path, String))
    cores = MN.decode_cores(data.cores)
    gauges = MN.decode_gauges(data.gauges)
    weights = haskey(data, :physical_weights) ? MN.parse_fraction.(data.physical_weights) : nothing
    if get(data, :source_adjoint, false)
        source_path = String(data.source_file)
        @assert bytes2hex(sha256(read(source_path))) == data.source_sha256
        source = JSON3.read(read(source_path, String))
        original = MN.decode_cores(source.candidate)
        @assert weights == MN.F[1, 4, 6, 4, 1]
        @assert cores == MN.adjoint_cores(original; weights)
        n = Int(source.length)
        @assert n == length(cores) == data.length && iseven(n)
        delta = MN.parse_fraction(source.beta)/(2big(2)^Int(source.power))
        radius = 5n*delta
        reference = (1+MN.parse_fraction(source.remainder)+radius+radius^2/2)/big(5)^div(n, 2)
        @assert MN.parse_fraction(data.analytic_reference_upper) == reference
    end
    bounds = MN.parse_fraction.(data.local_bounds)
    claimed = MN.parse_fraction(data.upper)
    @assert claimed == prod(bounds)
    @assert MN.verify_local_product(
        cores,
        gauges,
        bounds;
        precision = Int(data.precision_bits),
        weights,
    )
    println("Verified exact candidate MPO norm <= $(Float64(claimed)): $path")
    return claimed
end

if abspath(PROGRAM_FILE) == @__FILE__
    for path in ARGS
        verify_saved_norm(path)
    end
end
