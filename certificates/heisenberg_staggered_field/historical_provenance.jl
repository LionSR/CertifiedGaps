# Preserve the old certificate's meaning while validating the package rename.
# Included after common.jl. No historical JSON/hash is rewritten.
const HISTORICAL_HEAD = "1c4c07711072b8bfe1a60d188a498a491b8dc42b"
const MIGRATION_BASE = "22b2759e0790197b83e2571b465aa8dc69ebc0be"
function verify_migration_sources()
    count=0
    for h in FIELDS
        c=JSON3.read(read(joinpath(ROOT,"results","certificate_$(tag(h)).json"),String),Dict{String,Any})
        for (name,hash) in c["metadata"]["source_sha256"]
            archived=joinpath(ROOT,"historical_sources",name)
            @assert sha256_file(archived)==hash "historical source mismatch: $name"
            # All hashed executable sources differ only by the package identifier.
            current=startswith(name,"core/") ? joinpath(REPO,"julia",name[6:end]) : joinpath(ROOT,name)
            current=replace(current,"ThermalMPOGap"=>"CertifiedGaps")
            @assert read(current,String)==replace(read(archived,String),"ThermalMPOGap"=>"CertifiedGaps") "non-rename source change: $name"
            count+=1
        end
    end
    return count
end
const VERIFIED_MIGRATION_SOURCE_ENTRIES=verify_migration_sources()
