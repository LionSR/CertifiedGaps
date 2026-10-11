include("ladder.jl")
include("sparse_scan_basis.jl")
# Execute the identical scan driver with the sparse-to-real basis methods above.
include_string(Main,replace(read(joinpath(@__DIR__,"scan.jl"),String),
    "include(\"ladder.jl\")"=>""),joinpath(@__DIR__,"scan.jl"))
