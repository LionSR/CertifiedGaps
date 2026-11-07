include(joinpath(@__DIR__, "..", "src", "Certificate.jl"))
using .Certificate, LinearAlgebra
BLAS.set_num_threads(1)
g = isempty(ARGS) ? 24 : parse(Int, only(ARGS))
build_spectra(normpath(joinpath(@__DIR__, "..")), g)
