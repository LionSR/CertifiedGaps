# Run after installing the pinned toolchains. All files and logs stay in this directory.
using JSON, SHA
const HERE=@__DIR__
lake=get(ENV,"LAKE_EXE",joinpath(HERE,".tools/lean-4.35.0-rc3-linux/bin/lake"))
julia=joinpath(Sys.BINDIR,Base.julia_exename())
mkpath(joinpath(HERE,"logs"))
results=Dict{String,Any}()
cd(HERE) do
    for (name,cmd) in [("build",`$lake build`),
                       ("axioms",`$lake env lean Audit.lean`),
                       ("julia_tests",`$julia --project=. test/runtests.jl`)]
        log=joinpath(HERE,"logs","final_$name.log")
        seconds=@elapsed open(log,"w") do io
            run(pipeline(cmd;stdout=io,stderr=io))
        end
        results[name]=Dict("passed"=>true,"seconds"=>seconds)
    end
    audit=read(joinpath(HERE,"logs/final_axioms.log"),String)
    @assert !occursin("sorryAx",audit)
    @assert !occursin("ofReduceBool",audit)
    lines=filter(l->occursin("depends on axioms:",l),split(audit,'\n'))
    for l in lines
        a=match(r"\[(.*)\]",l)[1]
        @assert all(x->strip(x) in ["propext","Classical.choice","Quot.sound"],split(a,','))
    end
    @assert count(l->occursin("LeanSeedSmall.",l),split(audit,'\n'))==17
    # Check the protected source snapshot and exact pins independently of the build.
    manifest=JSON.parsefile("source_hashes.json")
    for (rel,h) in manifest["lean_sources"]
        original=read(joinpath(HERE,"../../lean",rel),String)
        copied=read(joinpath(HERE,"vendor",rel),String)
        @assert bytes2hex(sha256(original))==h
        if rel=="CertifiedGaps/Transfer/ProductLimit.lean"
            @assert split(copied,"/-!";limit=2)[2]==split(original,"/-!";limit=2)[2]
        else
            @assert copied==original
        end
    end
    cert=JSON.parsefile("certificate.json")
    @assert cert["source_hashes"]["certify.jl"]==bytes2hex(sha256(read("certify.jl")))
    @assert read("lean-toolchain",String)==read("../../lean/lean-toolchain",String)
    a=JSON.parsefile("lake-manifest.json")["packages"]
    b=JSON.parsefile("../../lean/lake-manifest.json")["packages"]
    @assert only(x["rev"] for x in a if x["name"]=="mathlib")==
            only(x["rev"] for x in b if x["name"]=="mathlib")
    results["axiom_audit"]="only the standard three; conditional hypotheses remain open"
    results["protected_sources"]="unchanged against recorded SHA256 snapshot"
    open("validation.json","w") do io; JSON.print(io,results,2); end
    println(JSON.json(results))
end
