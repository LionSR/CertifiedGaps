import CertifiedGaps
using CertifiedGaps.Results: write_result
using JSON3, SHA
const PG = CertifiedGaps.CertifiedMPO.PolynomialBounds

function polynomial_targets(
    output = joinpath(@__DIR__, "../../results/polynomial-thermal-targets.json"),
)
    rows = Any[]
    for power in (24, 26, 28, 30)
        c = PG.seed(256, 64, power)
        push!(
            rows,
            Dict(
                "length"=>c.n,
                "beta"=>string(c.beta),
                "squarings"=>power,
                "m"=>string(c.m),
                "delta"=>string(c.delta),
                "remainder_upper"=>string(c.remainder),
                "polynomial_coefficients"=>string.([c.identity, c.linear, c.quadratic]),
                "exact_polynomial_trace_inflation_upper"=>string(c.inflation_upper),
                "inflation_upper_display"=>Float64(c.inflation_upper),
                "scalar_majorant_certified"=>true,
                "thermal_trace_contracted"=>false,
            ),
        )
        println(
            "n=256 beta=64 squarings=$power exact-polynomial inflation <= $(Float64(c.inflation_upper))",
        )
    end
    result = Dict(
        "status"=>"scalar-polynomial-majorant-only",
        "complete"=>true,
        "uniform_gap_certified"=>false,
        "thermal_trace_certified"=>false,
        "targets"=>rows,
        "source_sha256"=>bytes2hex(
            sha256(read(joinpath(@__DIR__, "../../../../julia/src/CertifiedMPO/PolynomialBounds.jl"))),
        ),
    )
    write_result(output, result)
end

if abspath(PROGRAM_FILE) == @__FILE__
    polynomial_targets(
        isempty(ARGS) ? joinpath(@__DIR__, "../../results/polynomial-thermal-targets.json") :
        only(ARGS),
    )
end
