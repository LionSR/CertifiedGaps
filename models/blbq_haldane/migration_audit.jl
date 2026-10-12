using CertifiedGaps.Results, LinearAlgebra, TOML
BLAS.set_num_threads(1)
root=@__DIR__;repo=normpath(joinpath(root,"../.."))
oldref="ef0a3e6";rows=Any[]
for file in ("Results.jl","Seeds.jl","ExactDiagonalization.jl","Spins.jl",
    "ReflectionPositivity.jl","Recurrence.jl","recurrence/certified.jl",
    "recurrence/map.jl","recurrence/rational.jl")
    old=read(`git -C $repo show $(oldref*":julia/src/"*file)`,String)
    current=read(joinpath(repo,"julia/src",file),String)
    @assert replace(old,"ThermalMPOGap"=>"CertifiedGaps")==current
    push!(rows,(;file,identical_after_namespace_rename=true,
        current_sha256=sha256_file(joinpath(repo,"julia/src",file))))
end
oldmanifest=TOML.parse(read(`git -C $repo show $(oldref*":julia/Manifest.toml")`,String))
currentmanifest=TOML.parsefile(joinpath(repo,"julia/Manifest.toml"))
olddeps=oldmanifest["deps"];newdeps=currentmanifest["deps"]
for (name,entry) in olddeps
    name=="ThermalMPOGap" && continue
    @assert entry==newdeps[name] "dependency changed: $name"
end
write_result(joinpath(root,"results/migration_audit.json"),(;
    metadata=run_metadata((;old_core_ref=oldref,
        new_core_ref=strip(read(`git -C $repo rev-parse origin/main`,String)));
        source_sha256=source_hashes(root;dirs=(),files=("migration_audit.jl",)),
        status="core algorithms and dependency entries unchanged after namespace rename"),
    comparisons=rows,noncore_manifest_entries_identical=true,
    historical_evidence="historical/thermal_*.json and SHA-verified historical/thermal_source",
    note="Only certificate imports/project references migrated; no shared-core source edit."))
println("Migration audit passed: algorithms and dependency versions match")
