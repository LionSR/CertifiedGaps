# Run the Lean kernel checks; Julia never decides their truth.
include(joinpath(@__DIR__,"scripts","verify_hashes.jl"))
root=@__DIR__
lake=get(ENV,"LAKE_BIN",joinpath(root,".tools","elan","bin","lake"))
run(Cmd(`$lake build`;dir=root))
