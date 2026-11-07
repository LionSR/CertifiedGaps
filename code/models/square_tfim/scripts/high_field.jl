# This checks the arithmetic consequences of the graph theorem, without diagonalization.
include(joinpath(@__DIR__, "..", "src", "HighField.jl"))
using .HighField
using CertifiedGaps.ExactArithmetic: parse_rational

all(==("--check"), ARGS) || error("usage: high_field.jl [--check]")
root = dirname(@__DIR__)
"--check" in ARGS || certify_high_field(root)
record = check_high_field(root)
for row in record["points"]
    println("g/J=$(row["minimum_field"]): gap/J >= ",
        Float64(parse_rational(row["gap_lower"])))
end
