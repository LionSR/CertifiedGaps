include(joinpath(@__DIR__, "verify_operator_seed.jl"))
import ThermalMPOGap
const GS = ThermalMPOGap.CertifiedMPO.GaugeSearch
using LinearAlgebra
const MS = OE.MS
include(joinpath(@__DIR__, "operator_power_screen.jl"))
BLAS.set_num_threads(1)

"Compress an exact polynomial as A times the identity, preserving an operator-error proof."
function operator_seed_probe(
    source,
    output;
    maxdim = 16,
    bits = 64,
    candidate_arithmetic = :interval,
    conserve_charge = false,
)
    started = time()
    @assert abspath(source) != abspath(output)
    exact = checked_exact_seed(source)
    physical_charges = conserve_charge ? collect(-2:2) : nothing
    physical = sqrt.([Float64(PM.WEIGHTS[t]/PM.WEIGHTS[s]) for s in 1:5, t in 1:5])
    function gauges(cores)
        numerical = [Float64.(c) .* reshape(physical, 1, 5, 5, 1) for c in cores]
        proposed = GS.schmidt_gauge_proposal(numerical; physical_charges)
        return [OE.MS.relative_quantize(g; bits) for g in proposed]
    end
    first_gauges = gauges(exact.candidate)
    identity = PM.identity_mpo(exact.seed.n)
    println("Compressing exact polynomial at N=$(exact.seed.n)")
    flush(stdout)
    progress =
        row -> begin
            if row.site % 16 == 0 || row.site >= row.total-1
                println("Initial operator: $(row.stage), site $(row.site)/$(row.total)")
                flush(stdout)
            end
        end
    product = MP.product_candidate(
        exact.candidate,
        identity;
        maxdim,
        bits,
        weights = PM.WEIGHTS,
        first_gauges,
        candidate_arithmetic,
        physical_charges,
        progress,
    )
    @assert MP.verify_product(
        exact.candidate,
        identity,
        product;
        residual_arithmetic = :interval,
        progress,
    )
    println("Certifying initial operator error")
    flush(stdout)
    operator = OE.operator_error_certificate(
        exact.candidate,
        product;
        second_input = identity,
        first_gauges,
        prefix_gauges = gauges(product.cores),
        arithmetic = :interval,
        progress,
    )
    witness = splitext(output)[1]*"-witness.json"
    @assert abspath(witness) != abspath(source)
    open(witness, "w") do io
        JSON3.write(
            io,
            Dict(
                "product"=>MP.encode(exact.candidate, identity, product),
                "operator"=>OE.encode(operator),
            ),
        )
    end
    data = Dict(
        "construction"=>"operator-compressed-polynomial",
        "length"=>exact.seed.n,
        "source_file"=>source,
        "source_sha256"=>exact.hash,
        "witness_file"=>basename(witness),
        "witness_sha256"=>bytes2hex(sha256(read(witness))),
        "candidate"=>MN.encode_cores(product.cores),
        "operator_error_upper"=>string(operator.error),
        "hs_error_upper"=>string(product.error),
        "operator_error_display"=>Float64(operator.error),
        "hs_error_display"=>Float64(product.error),
        "candidate_arithmetic"=>string(candidate_arithmetic),
        "conserve_charge"=>conserve_charge,
        "generation_seconds"=>time()-started,
        "operator_error_over_reference_display"=>Float64(
            operator.error/PM.initial_operator_upper(exact.seed),
        ),
        "generator_replay_complete"=>false,
        "uniform_gap_certified"=>false,
        "thermal_trace_certified"=>false,
    )
    open(output, "w") do io
        JSON3.write(io, data)
    end
    verify_operator_seed(output)
    data["generator_replay_complete"] = true
    data["total_seconds"] = time()-started
    open(output, "w") do io
        JSON3.write(io, data)
    end
    return data
end

"Refine only the operator-error gauges of a saved linear seed compression."
function refine_operator_seed(source, output; powers = 0.0:0.05:0.6)
    @assert abspath(source) != abspath(output)
    source_text = read(source, String)
    data = JSON3.read(source_text, Dict{String,Any})
    witness_path = joinpath(dirname(source), data["witness_file"])
    witness_text = read(witness_path, String)
    @assert bytes2hex(sha256(codeunits(witness_text))) == data["witness_sha256"]
    raw = JSON3.read(witness_text, Dict{String,Any})
    A, identity, product = MP.decode(raw["product"])
    original = OE.decode(raw["operator"])
    screen = operator_power_screen(A, product, original; second_input = identity, powers)
    families = (original.first_gauges, original.second_gauges, original.prefix_gauges)
    changed = ntuple(k -> [s .* g for (s, g) in zip(screen.scales[k], families[k])], 3)
    operator = OE.operator_error_certificate(
        A,
        product;
        second_input = identity,
        first_gauges = changed[1],
        second_gauges = changed[2],
        prefix_gauges = changed[3],
        arithmetic = :interval,
    )
    raw["operator"] = OE.encode(operator)
    witness = splitext(output)[1]*"-witness.json"
    @assert abspath(witness) != abspath(witness_path)
    open(witness, "w") do io
        JSON3.write(io, raw)
    end
    data["witness_file"] = basename(witness)
    data["witness_sha256"] = bytes2hex(sha256(read(witness)))
    data["operator_error_upper"] = string(operator.error)
    data["operator_error_display"] = Float64(operator.error)
    exact = checked_exact_seed(data["source_file"])
    @assert exact.hash == data["source_sha256"]
    data["operator_error_over_reference_display"] =
        Float64(operator.error/PM.initial_operator_upper(exact.seed))
    data["refined_from"] = source
    data["refined_from_sha256"] = bytes2hex(sha256(codeunits(source_text)))
    data["selected_gauge_powers"] = screen.powers
    data["numerical_screen_certified"] = false
    data["screen_log_error_before"] = screen.baseline_log
    data["screen_log_error_after"] = screen.best_log
    data["generator_replay_complete"] = false
    open(output, "w") do io
        JSON3.write(io, data)
    end
    verify_operator_seed(output)
    data["generator_replay_complete"] = true
    open(output, "w") do io
        JSON3.write(io, data)
    end
    println(
        "Refined initial operator error: $(Float64(original.error)) -> $(Float64(operator.error))",
    )
    return data
end

if abspath(PROGRAM_FILE) == @__FILE__
    @assert length(ARGS) == 2
    operator_seed_probe(ARGS...)
end
