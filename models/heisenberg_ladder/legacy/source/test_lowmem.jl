include("ladder.jl")
include("sparse_scan_basis.jl")
include_string(Main,replace(read(joinpath(@__DIR__,"smoke.jl"),String),
    "include(\"ladder.jl\")"=>""),joinpath(@__DIR__,"smoke.jl"))
