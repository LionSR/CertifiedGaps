include(joinpath(@__DIR__, "src", "Certificate.jl"))
include(joinpath(@__DIR__, "src", "HighField.jl"))
using .Certificate, LinearAlgebra
using CertifiedGaps.ExactArithmetic: parse_rational
BLAS.set_num_threads(1)
args = filter(!=("--full"), ARGS)
points = isempty(args) ? POINTS : filter(p -> p.g == parse(Int, only(args)), POINTS)
isempty(points) && error("field is not a configured certificate point")
for p in points
    seconds = @elapsed cert = check(@__DIR__, p.g; full = "--full" in ARGS)
    println("Verified g=$(p.g): gap >= ",
        Float64(uniform_bound(cert)),
        " on every even square torus L >= 2; $seconds seconds")
end
HighField.check_high_field(@__DIR__)
println("Verified analytic Ising bound: gap >= 2sqrt(g(g-4J)) for every g > 4J")
