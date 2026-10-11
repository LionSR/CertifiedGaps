# Verify the rename and preserve all historical numerical and certified evidence.
using CertifiedGaps.Results
include("src/Provenance.jl")
using .Provenance
root=@__DIR__
repo=normpath(joinpath(root,"../.."))
snapshot=HISTORICAL_SNAPSHOT
paths=split(strip(read(`git -C $repo ls-tree -r --name-only $snapshot julia/src`,String)),'\n')
for path in paths
    old=read(`git -C $repo show $(snapshot*":"*path)`,String)
    current_path=replace(path,"ThermalMPOGap"=>"CertifiedGaps")
    @assert replace(old,"ThermalMPOGap"=>"CertifiedGaps")==read(joinpath(repo,current_path),String) path
end
old=read(`git -C $repo show $(snapshot*":certificates/two_leg_tfim/src/Ladder.jl")`,String)
@assert replace(old,"ThermalMPOGap"=>"CertifiedGaps")==read(joinpath(root,"src/Ladder.jl"),String)
result_paths=split(strip(read(`git -C $repo ls-tree -r --name-only $snapshot certificates/two_leg_tfim/results`,String)),'\n')
old_results=Dict{String,String}()
for path in result_paths
    endswith(path,".json") || continue
    relative=relpath(joinpath(repo,path),root)
    original=source_sha256(root,relative;snapshot)
    @assert original==source_sha256(root,relative) path
    old_results[relative]=original
end
write_result(joinpath(root,"results","migration_2026-10-11_provenance.json"),
    (;status="passed",historical_snapshot=snapshot,
      upstream_commit="22b2759e0790197b83e2571b465aa8dc69ebc0be",
      core_files_compared=length(paths),
      comparison="All core src bytes equal after ThermalMPOGap -> CertifiedGaps substitution and filename rename; Ladder.jl likewise.",
      preserved_historical_results=old_results,
      source_hashes=source_hashes(root;dirs=("src",),files=["check_migration.jl"])))
println("Namespace-only equality: $(length(paths)) core files and Ladder.jl; $(length(old_results)) historical JSON files unchanged.")
