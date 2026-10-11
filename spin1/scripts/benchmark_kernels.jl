import Spin1Certificates as SC
# Warm both methods before the timed same-machine comparison.
for kernel in (:reference, :optimized)
    SC.thermal_sector(1, 6, 0, 1; kernel)
    SC.thermal_sector(2, 6, 2, 1; kernel)
end
records = Any[]
for (source, n, twist, m) in ((1, 11, 1, 1), (2, 10, 2, 1))
    baseline = SC.thermal_sector(source, n, twist, m; kernel = :reference)
    optimized = SC.thermal_sector(source, n, twist, m; kernel = :optimized)
    @assert baseline.total == optimized.total
    push!(
        records,
        Dict(
            "source" => source,
            "n" => n,
            "twist" => twist,
            "magnetization" => m,
            "reference_seconds" => baseline.seconds,
            "optimized_seconds" => optimized.seconds,
            "speedup" => baseline.seconds / optimized.seconds,
            "total" => string(baseline.total),
        ),
    )
    println(
        "kernel source=$source n=$n reference=$(baseline.seconds)s optimized=$(optimized.seconds)s speedup=$(baseline.seconds/optimized.seconds)",
    )
    flush(stdout)
end
SC.write_result(
    joinpath(@__DIR__, "../results/kernel-benchmark.json"),
    Dict("status" => "passed", "comparisons" => records),
)
