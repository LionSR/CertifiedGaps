# Recombine the existing source-1 moments and length-120 trial certificate.
# Discovery may use floating point; every accepting inequality below is rational.
include(joinpath(@__DIR__, "optimize_bounds.jl"))

using CertifiedGaps.ExactArithmetic: integer_root, power_interval, power_upper

"Source-1 upper bounds at every real moment order k >= 12."
function source1_moment_bounds(k::Int, inputs)
    @assert k >= 12
    U, V, mJ, mN = (inputs[key] for key in ("U", "V", "mJ", "mN"))
    @assert U^12 < mJ < 2U^12 && V^12 < mN < 2V^12
    j = U^k + power_upper(mJ - U^12, frac(k, 12))
    n = V^k + power_upper(mN - V^12, frac(k, 12))
    return j + 3n, j
end

"Upper bound on the residual s-moment after removing a J-sector weight > t."
function lifted_residual(D::RationalBig, tau::RationalBig, t::RationalBig, s::RationalBig)
    @assert 0 < t <= tau <= D && s >= 1
    at_zero = 3power_upper((D - t) / 3, s)
    at_maximum = power_upper(tau - t, s) + 3power_upper((D - tau) / 3, s)
    return max(at_zero, at_maximum)
end

function lifted_seed(inputs; k::Int, multiplier::RationalBig, trial_length::Int = 120)
    @assert iseven(k) && iseven(trial_length) && 12 <= k <= trial_length / 2
    @assert multiplier >= 1
    moment, twisted = source1_moment_bounds(k, inputs)
    D = power_upper(moment, multiplier)
    E = power_upper(twisted, multiplier)
    tau = (D + 3E) / 4
    h = frac(trial_length, 2k)
    exponent = 2h
    # If all weights were <= t, their exponent-moment would be <= D*t^(exponent-1).
    # A weight > t must be in J because an N weight occurs three times.
    t = first(power_interval(inv(D), inv(exponent - 1))) - frac(1, big(10)^24)
    @assert 0 < t < tau && 3t > D
    @assert D * power_upper(t, exponent - 1) < 1
    initial_t = t
    refinements = Any[]
    for _ in 1:40
        residual = lifted_residual(D, tau, t, exponent)
        @assert residual < 1
        s = first(power_interval(1 - residual, inv(exponent))) - frac(1, big(10)^24)
        s <= t && break
        @assert s < tau && power_upper(s, exponent) + residual < 1
        push!(
            refinements,
            Dict(
                "old_lower" => string(t),
                "new_lower" => string(s),
                "residual" => string(residual),
            ),
        )
        t = s
    end
    residual = lifted_residual(D, tau, t, h)
    dominant = first(power_interval(t, h))
    p = 1 - (1 + residual / dominant)^(-2)
    total, _ = source1_moment_bounds(trial_length, inputs)
    @assert total > 1
    excited_ratio = power_upper(total - 1, multiplier)
    q = 1 - (1 + excited_ratio)^(-2)
    details = Dict(
        "moment_length" => k,
        "trial_length" => trial_length,
        "temperature_multiplier" => string(multiplier),
        "moment_upper" => string(moment),
        "twisted_moment_upper" => string(twisted),
        "D" => string(D),
        "E" => string(E),
        "tau" => string(tau),
        "h" => string(h),
        "initial_dominant_lower" => string(initial_t),
        "dominant_lower" => string(t),
        "refinements" => refinements,
        "spatial_residual_upper" => string(residual),
        "spatial_dominant_power_lower" => string(dominant),
        "physical_partition_upper" => string(total),
        "physical_excited_ratio_upper" => string(excited_ratio),
        "p" => string(p),
        "q" => string(q),
    )
    return p, q, details
end

function lift_main(output = joinpath(@__DIR__, "temperature-lift.json"))
    SC.write_result(output, Dict("status" => "running"))
    started = time()
    provenance = SC.verify_source()
    reproduction_path = joinpath(@__DIR__, "../results/verification.json")
    reproduction = SC.read_exact_json(reproduction_path)
    @assert reproduction["status"] == "passed"
    @assert reproduction["implementation_sha256"] == SC.implementation_fingerprint()
    for (name, component) in reproduction["components"]
        @assert bytes2hex(sha256(read(joinpath(@__DIR__, "../results", "$name.json")))) ==
                component["sha256"]
    end
    oldp, oldq, input_strings = tightened_source1(SC.source1_initialization())
    inputs = Dict(key => SC.parse_rational(input_strings[key]) for key in ("U", "V", "mJ", "mN"))
    # Independently check the two finite variational signs used below.
    trial_hashes = Dict{String,String}()
    for n in (72, 120)
        path = joinpath(SC.SOURCE, "evidence/trials/certificate$n.json")
        c = SC.read_exact_json(path)
        U = SC.integer(c[n == 72 ? "one_bond_energy_numerator" : "U"])
        V = SC.integer(c[n == 72 ? "norm" : "V"])
        @assert V > 0 && 500000U + 700741V < 0
        trial_hashes[string(n)] = bytes2hex(sha256(read(path)))
    end
    A = (1 - oldp)^2 * (1 - oldq)
    bridge_lower = first(log_enclosure(A / (1 - A))) / frac(49, 4)
    records = Dict{String,Any}()
    for (name, k, multiplier) in (("simple", 30, frac(2)), ("optimized", 28, frac(17, 8)))
        p, q, details = lifted_seed(inputs; k, multiplier)
        # The stronger seeds square below the older 10^-512 grid too early.
        # Keep the rounding floor well below both defects through the tail step.
        chain =
            rational_gap_chain(p, q, 60, frac(49, 4) * multiplier; steps = 10, grid = big(10)^2048)
        claim = SC.parse_rational(chain["gap_strict_lower"])
        @assert 0 < claim < bridge_lower
        records[name] = Dict(
            "initialization" => details,
            "chain" => chain,
            "uniform_even_length_threshold" => 36,
            "uniform_gap_strict_lower" => string(claim),
            "uniform_gap_strict_lower_decimal" => Float64(claim),
        )
        println(
            "$name: gamma > $(Float64(claim)), every even L >= 36; p=$(Float64(p)), q=$(Float64(q))",
        )
        flush(stdout)
    end
    SC.write_result(
        output,
        Dict(
            "status" => "exact-scalar-checks-passed",
            "results" => records,
            "inputs" => input_strings,
            "bridge" => Dict(
                "length_from" => 36,
                "length_to" => 72,
                "beta" => "49//4",
                "p" => string(oldp),
                "q" => string(oldq),
                "purity_lower" => string(A),
                "gap_lower" => string(bridge_lower),
            ),
            "source_integrity" => provenance,
            "reproduction_verification_sha256" => bytes2hex(sha256(read(reproduction_path))),
            "trial_sha256" => trial_hashes,
            "script_sha256" => bytes2hex(sha256(read(@__FILE__))),
            "bootstrap_script_sha256" =>
                bytes2hex(sha256(read(joinpath(@__DIR__, "optimize_bounds.jl")))),
            "new_thermal_or_trial_contractions" => false,
            "execution_mode" => "fresh-rational-temperature-lifting-on-reproduced-finite-inputs",
            "scope" => "Exact finite scalar arithmetic plus the documented temperature-lifting lemma and original analytic spatial-transfer and bootstrap theorems; not a formal proof.",
            "julia" => string(VERSION),
            "completed_utc" => string(now(UTC)),
            "wall_seconds" => time() - started,
        ),
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    lift_main(isempty(ARGS) ? joinpath(@__DIR__, "temperature-lift.json") : only(ARGS))
end
