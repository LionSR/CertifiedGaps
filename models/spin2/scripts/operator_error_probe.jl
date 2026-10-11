import CertifiedGaps
using CertifiedGaps.Results: write_result
using JSON3, SHA, LinearAlgebra
using CertifiedGaps.CertifiedMPO.GaugeSearch
const OE = CertifiedGaps.CertifiedMPO.OperatorError
const GS = CertifiedGaps.CertifiedMPO.GaugeSearch
const MS = OE.MS
const MN = OE.MN
BLAS.set_num_threads(1)

function operator_probe(
    source,
    output;
    iterations = 200,
    gauge_mode = :diagonal,
    gauge_powers = (0.5, 0.5, 0.5),
    arithmetic = :interval,
    precision = 192,
)
    @assert abspath(source) != abspath(output)
    @assert gauge_mode in (:diagonal, :schmidt) && length(gauge_powers) == 3
    source_text = read(source, String)
    data = JSON3.read(source_text)
    input, square = MS.decode_witness(data)
    weights = square.weights
    d = size(first(input), 2)
    physical =
        isnothing(weights) ? ones(d, d) :
        sqrt.([Float64(weights[t]/weights[s]) for s in 1:d, t in 1:d])
    function search(cores, power)
        numerical = [Float64.(c) .* reshape(physical, 1, d, d, 1) for c in cores]
        if gauge_mode == :schmidt
            proposed = GS.schmidt_gauge_proposal(numerical)
            scaled = [g .* vec(sqrt.(sum(abs2, g; dims = 2))) .^ (2power-1) for g in proposed]
            println("Schmidt operator-gauge proposal, power=$power")
            flush(stdout)
            return [MS.relative_quantize(g; bits = 64) for g in scaled]
        end
        @assert power == 0.5
        found = optimize_gauge(numerical; iterations)
        println("Gauge product estimate: $(exp(found.initial)) -> $(exp(found.final))")
        flush(stdout)
        return MN.rational_gauges(cores, found.x)
    end
    started = time()
    first_copy = square.adjoint_first ? MN.adjoint_cores(input; weights) : input
    progress = row -> begin
        if row.site % 16 == 0 || row.site == row.total-1
            println("Operator residual: site $(row.site)/$(row.total-1)")
            flush(stdout)
        end
    end
    proof = OE.operator_error_certificate(
        input,
        square;
        first_gauges = search(first_copy, gauge_powers[1]),
        second_gauges = search(input, gauge_powers[2]),
        prefix_gauges = search(square.cores, gauge_powers[3]),
        arithmetic,
        precision,
        progress,
    )
    encoded = OE.encode(proof)
    encoded["source_file"] = source
    encoded["source_sha256"] = bytes2hex(sha256(codeunits(source_text)))
    encoded["operator_error_upper_display"] = Float64(proof.error)
    encoded["hs_error_upper_display"] = Float64(square.error)
    encoded["seconds"] = time()-started
    encoded["uniform_gap_certified"] = false
    encoded["thermal_trace_certified"] = false
    encoded["gauge_mode"] = string(gauge_mode)
    encoded["gauge_powers"] = collect(gauge_powers)
    encoded["precision_bits"] = precision
    encoded["generator_replay_complete"] = false
    write_result(output, encoded)
    println("Replaying operator residual certificate")
    flush(stdout)
    @assert OE.verify_operator_error(
        input,
        square,
        proof;
        residual_arithmetic = arithmetic,
        precision,
        progress,
    )
    encoded["generator_replay_complete"] = true
    encoded["seconds"] = time()-started
    write_result(output, encoded)
    println(
        "Operator error <= $(Float64(proof.error)); HS bound=$(Float64(square.error)); seconds=$(time()-started)",
    )
    return encoded
end

if abspath(PROGRAM_FILE) == @__FILE__
    operator_probe(ARGS...)
end
