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
    println("g=$(p.g) beta=$beta defects=$(cert["seed"]["defects"])")
    println("Uniform square-torus gap >= ",
        Float64(parse_rational(cert["chain"]["gap_lower_bound"])))
    println("Elapsed: $seconds seconds")
end
