import CertifiedGaps
using CertifiedGaps.Results: write_result
using JSON3, SHA
import Nemo
const PG = CertifiedGaps.CertifiedMPO.PolynomialBounds

"Compare rigorously bounded polynomial overhead, without an MPO claim."
function degree_screen(
    output;
    n = 256,
    beta = 64,
    degrees = (2, 4, 6, 8, 10, 12, 16, 20, 24, 32),
    max_inflation = 1001//1000,
    max_power = 32,
)
    RR = Nemo.ArbField(384)
    log_budget = log(RR(max_inflation))
    rows = Any[]
    for degree in degrees
        selected = nothing
        for power in 0:max_power
            seed = PG.taylor_seed(n, beta, power, degree)
            # Compare logarithms: large rejected values need not be exponentiated.
            if RR(seed.log_inflation_upper) < log_budget
                selected = seed
                break
            end
        end
        @assert !isnothing(selected)
        c = selected
        x = RR(c.delta*c.spectral_radius)
        @assert RR(c.remainder) > exp(x)*x^(degree+1)/factorial(big(degree+1))
        @assert RR(c.log_inflation_upper) > 4c.m*RR(c.remainder)*exp(x)
        inflation = PG.relative_upper(exp(RR(c.log_inflation_upper)))
        @assert inflation < max_inflation
        row = Dict(
            "degree"=>degree,
            "power"=>c.power,
            "delta"=>string(c.delta),
            "scaled_spectral_radius"=>string(c.delta*c.spectral_radius),
            "remainder_upper"=>string(c.remainder),
            "log_inflation_upper"=>string(c.log_inflation_upper),
            "inflation_upper"=>string(inflation),
            "inflation_upper_display"=>Float64(inflation),
            "coefficients"=>string.(c.coefficients),
            "polynomial_operator_upper"=>string(c.operator_upper),
        )
        push!(rows, row)
        println("degree=$degree squarings=$(c.power) inflation<=$(Float64(inflation))")
        flush(stdout)
    end
    result = Dict(
        "status"=>"certified-scalar-polynomial-overhead-only",
        "scalar_majorant_certified"=>true,
        "compressed_mpo_certified"=>false,
        "thermal_trace_certified"=>false,
        "uniform_gap_certified"=>false,
        "length"=>n,
        "beta"=>string(beta),
        "max_inflation"=>string(max_inflation),
        "verification_precision_bits"=>384,
        "source_sha256"=>bytes2hex(sha256(read(@__FILE__))),
        "bounds_source_sha256"=>bytes2hex(
            sha256(read(joinpath(@__DIR__, "../../../julia/src/CertifiedMPO/PolynomialBounds.jl"))),
        ),
        "rows"=>rows,
    )
    write_result(output, result)
    return result
end

if abspath(PROGRAM_FILE) == @__FILE__
    @assert length(ARGS) == 1
    degree_screen(only(ARGS))
end
