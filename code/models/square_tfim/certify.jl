include(joinpath(@__DIR__, "src", "Certificate.jl"))
using .Certificate, LinearAlgebra
using CertifiedGaps.ExactArithmetic: parse_rational
BLAS.set_num_threads(1)
args = filter(!=("--reuse"), ARGS)
points = isempty(args) ? POINTS : filter(p -> p.g == parse(Int, args[1]), POINTS)
isempty(points) && error("field is not a configured certificate point")
length(args) <= 2 || error("usage: certify.jl [g [beta]] [--reuse]")
for p in points
    beta = length(args) < 2 ? p.beta : parse_rational(args[2])
    seconds = @elapsed cert = certify(@__DIR__, p.g, beta; reuse = "--reuse" in ARGS)
    direct = cert["direct_temperature"]
    println("g=$(p.g), direct inputs: ", direct["input"]["defects"])
    println("Side-four seed beta=$(direct["seed"]["beta"]), defects=",
        direct["seed"]["defects"])
    println("Uniform square-torus gap >= ", Float64(uniform_bound(cert)), " for even L >= 4")
    region = cert["parameter_interval"]
    println("Field interval ", region["field_interval"], ": gap >= ",
        Float64(parse_rational(region["gap_lower_bound"])))
    println("Elapsed: $seconds seconds")
end
