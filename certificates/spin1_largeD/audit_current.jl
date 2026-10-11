# Extra audit checks for replaying historical certificates with the renamed current core.
using SHA
const historical_core_revision = "ef0a3e6cf54ec56dcef5e8a37963d0184f6fd31b"

function recorded_source_matches(file,hash)
    any(base->isfile(joinpath(base,file)) && sha256_file(joinpath(base,file))==hash,
        (root,joinpath(root,"provenance","legacy")))
end

function recorded_core_matches(file,hash)
    current=joinpath(root,"../../julia",file)
    isfile(current) && sha256_file(current)==hash && return true
    repository=normpath(joinpath(root,"../.."))
    object=historical_core_revision*":julia/"*file
    bytes=read(`git -C $repository show $object`)
    bytes2hex(sha256(bytes))==hash
end

function rename_only(file)
    old=read(joinpath(root,"provenance","legacy",file),String)
    replace(old,"ThermalMPOGap"=>"CertifiedGaps")==read(joinpath(root,file),String)
end
