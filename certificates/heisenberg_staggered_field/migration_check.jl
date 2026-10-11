# Validate the separately recorded migration replay and preserved old evidence.
include("common.jl")
include("historical_provenance.jl")
old=JSON3.read(read(joinpath(ROOT,"results","spectral_replay.json"),String),Dict{String,Any})
freshpath=joinpath(ROOT,"results","migration_spectral_replay.json")
fresh=JSON3.read(read(freshpath,String),Dict{String,Any})
@assert fresh["status"]=="complete" && length(fresh["records"])==5
accepted=Any[]
for a in old["records"]
    b=only(filter(r->r["h"]==a["h"],fresh["records"]))
    s=enclosed_seeds(b["small_sectors"],b["large_sectors"],rat(b["beta"]))
    @assert s.s==rat(b["epsilon_s"]) && s.t==rat(b["epsilon_t"])
    gap=certified_gap_bound(s.s,s.t,rat(b["beta"])).bound
    @assert gap==b["bound"]
    published=rat(a["published_bound_rational"])
    @assert Rational{BigInt}(gap)>=published>0
    push!(accepted,Dict("h"=>a["h"],"bound"=>string(published),"fresh_raw_bound"=>gap))
end
manifest=JSON3.read(read(joinpath(ROOT,"results","audit_manifest.json"),String),Dict{String,Any})
preserved=0
for (path,hash) in manifest["file_sha256"]
    startswith(path,"results/") || continue
    @assert sha256_file(joinpath(ROOT,path))==hash
    global preserved+=1
end
write_result(joinpath(ROOT,"results","migration_check.json"),Dict(
    "status"=>"passed","accepted"=>accepted,"historical_head"=>HISTORICAL_HEAD,
    "migration_base"=>MIGRATION_BASE,"historical_results_verified"=>preserved,
    "source_entries_verified"=>VERIFIED_MIGRATION_SOURCE_ENTRIES,
    "migration_replay_sha256"=>sha256_file(freshpath),
    "source_sha256"=>file_hashes(ROOT,["common.jl","historical_provenance.jl","migration_check.jl"])))
println("All five conservative bounds supported; ",preserved,
    " historical result hashes and ",VERIFIED_MIGRATION_SOURCE_ENTRIES," source entries verified.")
