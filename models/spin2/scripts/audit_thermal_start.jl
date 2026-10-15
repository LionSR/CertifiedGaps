using Spin2Gap
import CertifiedGaps
import Nemo
using CertifiedGaps.Results: source_hashes, run_metadata, write_result

const ThermalStartBounds = CertifiedGaps.CertifiedMPO.PolynomialBounds

"Return strict rational endpoints of an Arb ball, independently checked in its parent field."
function audit_endpoints(value)
    upper = ThermalStartBounds.relative_upper(value)
    lower = -ThermalStartBounds.relative_upper(-value)
    field = parent(value)
    @assert field(lower) < value < field(upper)
    return Dict("lower" => string(lower), "upper" => string(upper))
end

"""
    audit_thermal_start()

Replay the scalar starting-scale and four-site Trotter tables in `IDEAS.md` using Arb.
Checks cover polynomial remainders, the exact four-site spectrum with all multiplicities,
and uncompressed Trotter traces. They do not certify a compressed MPO or a large-ring seed.
"""
function audit_thermal_start()
    started = time_ns()
    field = Nemo.ArbField(384)
    budget = 1001//1000
    degrees = (2, 4, 8, 12, 16, 32, 64, 128, 256, 512)
    expected_powers = (28, 21, 17, 16, 15, 13, 12, 11, 10, 9)
    polynomial_rows = Any[]
    for (degree, expected_power) in zip(degrees, expected_powers)
        selected = nothing
        for power in 0:32
            seed = ThermalStartBounds.taylor_seed(256, 64, power, degree)
            # Compare logarithms so rejected large bounds never need exponentiation.
            if field(seed.log_inflation_upper) < log(field(budget))
                selected = seed
                break
            end
        end
        @assert !isnothing(selected)
        seed = selected
        @assert seed.power == expected_power
        x = field(seed.delta*seed.spectral_radius)
        @assert field(seed.remainder) > exp(x)*x^(degree+1)/factorial(big(degree+1))
        @assert field(seed.log_inflation_upper) > 4seed.m*field(seed.remainder)*exp(x)
        inflation = ThermalStartBounds.relative_upper(exp(field(seed.log_inflation_upper)))
        @assert inflation < budget
        push!(
            polynomial_rows,
            Dict(
                "degree" => string(degree),
                "squarings" => string(seed.power),
                "starting_beta" => string(2seed.delta),
                "remainder_upper" => string(seed.remainder),
                "log_inflation_upper" => string(seed.log_inflation_upper),
                "inflation_upper" => string(inflation),
            ),
        )
        println("degree=$degree squarings=$(seed.power) starting_beta=$(2seed.delta)")
    end

    field = Nemo.ArbField(192)
    beta, shift = 4, 238061//50000
    gibbs, dimension = field(0), 0
    # H4 = (S1+S3) dot (S2+S4); couple the two diagonal-pair spins u and v.
    for u in 0:4, v in 0:4, s in abs(u-v):(u+v)
        numerator = s*(s+1)-u*(u+1)-v*(v+1)
        @assert iseven(numerator)
        energy, multiplicity = div(numerator, 2), 2s+1
        gibbs += multiplicity*exp(-field(beta)*(energy+4field(shift)))
        dimension += multiplicity
    end
    @assert dimension == 5^4
    @assert gibbs > 0
    trotter_rows = Any[]
    previous_ratio = nothing
    for steps in (64, 128, 256)
        reference = Spin2Gap.ThermalBounds.trotter_upper(4, beta, steps; shift)
        counted = sum((s.magnetization == 0 ? 1 : 2)*s.dimension for s in reference.sectors)
        @assert counted == dimension
        ratio = reference.trotter_trace/gibbs
        @assert ratio > 1
        isnothing(previous_ratio) || @assert ratio < previous_ratio
        previous_ratio = ratio
        push!(
            trotter_rows,
            Dict(
                "steps" => string(steps),
                "local_step" => string(beta//steps),
                "trotter_trace" => audit_endpoints(reference.trotter_trace),
                "trotter_over_gibbs" => audit_endpoints(ratio),
            ),
        )
        println("steps=$steps Trotter/Gibbs=$ratio")
        flush(stdout)
    end

    root = normpath(joinpath(@__DIR__, "../../.."))
    hashes = source_hashes(
        root;
        dirs = ("julia/src", "models/spin2/src"),
        files = (
            relpath(@__FILE__, root),
            "models/spin2/Project.toml",
            "models/spin2/Manifest.toml",
        ),
    )
    config = Dict(
        "polynomial_length" => "256",
        "polynomial_beta" => "64//1",
        "max_inflation" => string(budget),
        "polynomial_check_bits" => "384",
        "trotter_check_bits" => "192",
        "reference_length" => "4",
        "reference_beta" => string(beta//1),
        "reference_shift" => string(shift),
    )
    return Dict(
        "status" => "passed-finite-thermal-start-audit",
        "uniform_gap_certified" => false,
        "large_ring_thermal_seed_certified" => false,
        "compressed_mpo_certified" => false,
        "source_hashes" => hashes,
        "run_metadata" => run_metadata(
            config;
            source_sha256 = hashes[relpath(@__FILE__, root)],
            packages = (Spin2Gap, CertifiedGaps, Nemo),
            status = "finite-diagnostics-only",
            certified = false,
        ),
        "polynomial_rows" => polynomial_rows,
        "four_site_gibbs_trace" => audit_endpoints(gibbs),
        "trotter_rows" => trotter_rows,
        "elapsed_nanoseconds" => string(time_ns()-started),
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) <= 1 || error("usage: audit_thermal_start.jl [output.json]")
    result = audit_thermal_start()
    isempty(ARGS) || write_result(only(ARGS), result)
    println("Audit passed; no large-ring thermal seed or uniform gap is certified.")
end
