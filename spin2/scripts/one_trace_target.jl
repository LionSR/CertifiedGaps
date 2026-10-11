# Conditional scalar targets only. No spin-2 thermal upper bound is assumed to
# have been computed by this script.
using ThermalMPOGap.Results: write_result
module OneTraceTargets
using ThermalMPOGap.ExactArithmetic: RationalBig, frac, first_true, parse_rational
using ThermalMPOGap.Recurrence: rational_gap_chain
const F = RationalBig
# The exact scalar chain was spin1's `experiments/optimize_bounds.jl` `gap_chain`; it is now the
# verbatim core port (default grid 10^512, same keywords and output keys).
const gap_chain = rational_gap_chain
const SCALAR_BACKEND = joinpath(@__DIR__, "../../julia/src/recurrence/rational.jl")

function target(U::F, n::Int, beta::F; h = 2, r = 3)
    @assert U > 1 && iseven(n) && n >= 4 && h >= 1 && r >= 1
    D = U^r
    D < 2 || return nothing
    t, _ = first_true(x->D*x^(2h-1) >= 1, frac(0), frac(1))
    D*t^(2h-1) < 1 && 2t > D || return nothing
    if t^(2h)+(D-t)^(2h) < 1
        t, _ = first_true(x->x^(2h)+(D-x)^(2h) >= 1, t, frac(1))
    end
    @assert t^(2h)+(D-t)^(2h) < 1
    p = 1-(1+((D-t)/t)^h)^(-2)
    Q = U^(2h)
    q = 1-(1+(Q-1)^r)^(-2)
    proof = try
        gap_chain(p, q, h*n, r*beta; steps = 8)
    catch error
        error isa AssertionError || rethrow()
        return nothing
    end
    return Dict(
        "status"=>"conditional-scalar-target-only",
        "thermal_length"=>n,
        "thermal_beta"=>string(beta),
        "thermal_shifted_trace_upper_target"=>string(U),
        "trial_length_required"=>2h*n,
        "temperature_multiplier"=>r,
        "spatial_multiplier"=>h,
        "dominant_spatial_weight_lower"=>string(t),
        "spatial_purity_lower"=>string(1-p),
        "temporal_purity_lower"=>string(1-q),
        "conditional_gap_proof"=>proof,
        "spin2_thermal_input_certified"=>false,
    )
end
end

using .OneTraceTargets, JSON3, SHA
const OT = OneTraceTargets
function main(
    output = joinpath(@__DIR__, "../results/one-trace-targets.json");
    witness = joinpath(@__DIR__, "../results/trial-137-certificate.json"),
    n = 128,
    beta = OT.frac(32),
    shift = OT.frac(238061, 50000),
)
    trial = JSON3.read(read(witness, String), Dict{String,Any})
    @assert trial["complete"] && trial["status"] == "finite-energy-certificates-passed"
    upper(c) =
        parse(BigInt, c["energy_per_bond_interval"]["upper_numerator"])//parse(
            BigInt,
            c["energy_per_bond_interval"]["denominator"],
        )
    finite = only(filter(c->c["length"]==n, trial["certificates"]))["energy_per_bond_interval"]
    trial_upper = parse(BigInt, finite["upper_numerator"])//parse(BigInt, finite["denominator"])
    exponent_lower = -beta*n*(shift+trial_upper)
    @assert exponent_lower > 0
    records = Any[]
    for num in (102, 104, 106, 108, 110, 112)
        U = OT.frac(num, 100)
        candidates = Any[]
        for h in 1:4, r in 1:5
            any(t->t["length"]==2h*n && upper(t) < -shift, trial["certificates"]) || continue
            c = OT.target(U, n, beta; h, r)
            isnothing(c) || push!(candidates, c)
        end
        if isempty(candidates)
            println("U=$(Float64(U)): no passing target in scanned range")
        else
            best = candidates[argmax([
                c["conditional_gap_proof"]["gap_strict_lower_decimal"] for c in candidates
            ])]
            # e^x > 1+x for x>0 is an exact rational rejection screen.
            best["current_certified_shift"] = string(shift)
            best["required_trial_condition_certified"] = true
            best["variational_trace_lower_strict"] = string(1+exponent_lower)
            best["thermal_target_ruled_out_at_current_shift"] = 1+exponent_lower >= U
            push!(records, best)
            println(
                "U=$(Float64(U)): h=$(best["spatial_multiplier"]) r=$(best["temperature_multiplier"]) gap=$(best["conditional_gap_proof"]["gap_strict_lower_decimal"]) trial L=$(best["trial_length_required"])",
            )
        end
    end
    write_result(
        output,
        Dict(
            "status"=>"conditional-targets",
            "uniform_gap_certified"=>false,
            "targets"=>records,
            "source_sha256"=>bytes2hex(sha256(read(@__FILE__))),
            "scalar_backend_sha256"=>bytes2hex(sha256(read(OT.SCALAR_BACKEND))),
            "variational_screen_source_sha256"=>bytes2hex(sha256(read(witness))),
        ),
    )
end
if abspath(PROGRAM_FILE) == @__FILE__
    main(isempty(ARGS) ? joinpath(@__DIR__, "../results/one-trace-targets.json") : only(ARGS))
end
