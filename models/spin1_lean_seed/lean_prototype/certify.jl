# Generate arithmetic witness recipes. This is NOT a complete gap certificate.
include(joinpath(@__DIR__,"scripts","generate.jl"))
for (n,q) in ((6,24),(8,32),(10,16))
    generate(n,q)
end
println("Generated untrusted recipes. Run check.jl with the pinned Lean environment.")
