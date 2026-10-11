# One-time, provenance-preserving migration. Run only after both legacy computations finish.
if Base.find_package("ThermalMPOGap")!==nothing
    @eval using ThermalMPOGap.Results
else
    @eval using CertifiedGaps.Results
end
root=@__DIR__
for suffix in ("4_1","1603_500")
    @assert isfile(joinpath(root,"results","certificate_D$suffix.json"))
end
legacy=joinpath(root,"provenance","legacy")
@assert !isdir(legacy) "legacy sources have already been archived"
files=["model.jl","certify.jl","certify_parallel.jl","scan.jl","refine_beta.jl",
    "seed.jl","report.jl","survey_check.jl","check.jl","Project.toml",
    "test/runtests.jl","test/small.jl","test/symmetry.jl"]
changes=Dict{String,Any}()
for file in files
    old=joinpath(root,file); saved=joinpath(legacy,file)
    mkpath(dirname(saved)); cp(old,saved)
    original=read(old,String)
    migrated=replace(original,"ThermalMPOGap"=>"CertifiedGaps")
    write(old,migrated)
    changes[file]=Dict("historical_sha256"=>sha256_file(saved),
        "renamed_sha256"=>sha256_file(old))
end
write_result(joinpath(root,"results","package-migration.json"),Dict(
    "historical_core"=>"ef0a3e6cf54ec56dcef5e8a37963d0184f6fd31b",
    "target_core"=>"22b2759e0790197b83e2571b465aa8dc69ebc0be",
    "transformation"=>"literal package identifier replacement only; checker integration is validated separately",
    "files"=>changes,"source_hashes"=>source_hashes(root;dirs=(),files=("migrate_package.jl",))))
println("Archived historical sources; migrated package identifiers inside the certificate only.")
