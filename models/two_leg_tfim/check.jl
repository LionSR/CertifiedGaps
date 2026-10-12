# Fast checker by default. --full repeats the expensive independently reproducible eigensolve.
if "--full" in ARGS
    filter!(!=("--full"),ARGS)
    isempty(ARGS) && push!(ARGS,"3//1")
    include("certify.jl")
end
include("test/runtests.jl")
