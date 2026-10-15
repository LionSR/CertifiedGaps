include("common.jl")
using Downloads
using Nemo
const SOURCE=joinpath(@__DIR__,"upstream","verification/computations")
integer(x::AbstractString)=parse(BigInt,x)
integer(x::Integer)=BigInt(x)
include(joinpath(ROOT,"models/spin1/src/trials.jl"))
function fetch_matrix()
    base="https://raw.githubusercontent.com/openai/math/" *
        "adc7f1241b42e322a6451854ab7e4b4c146bf78a/preprints/" *
        "The-periodic-spin-one-Haldane-gap-September-24-2026/verification/computations"
    mkpath(joinpath(SOURCE,"evidence/trials"))
    man=joinpath(SOURCE,"source-manifest.json")
    isfile(man) || Downloads.download("$base/source-manifest.json",man)
    @assert sha256_file(man)=="2d2461e8f87ef9c3af62ce226412b106fcfb34f01d5139683a90708859e2dc92"
    manifest=read_exact_json(man)
    key="evidence/trials/matrix-data.json"
    expected=only(v[key] for (k,v) in manifest if v isa AbstractDict && haskey(v,key))
    path=joinpath(SOURCE,key)
    isfile(path) || Downloads.download("$base/$key",path)
    @assert sha256_file(path)==expected
end
function trial_main()
    fetch_matrix()
    labels,A=trial_matrices()
    blocks=trial_blocks(labels,A;native=true)
    physical_trial_check(A,blocks)
    scale=maximum(maximum(abs,eigvals(Float64.(b.T))) for b in blocks)
    floating=Any[]
    for n in 24:4:192
        u=v=0.0
        for b in blocks
            T=Float64.(b.T)/scale;O=Float64.(b.O)/scale^2
            X=T^(n-2)
            u+=sum(X.*transpose(O));v+=sum(X.*transpose(T*T))
        end
        push!(floating,Dict("length"=>n,"energy_per_bond"=>u/v,"certified"=>false))
    end
    started=time()
    n=isempty(ARGS) ? 144 : parse(Int,first(ARGS))
    u=v=big(0)
    records=Any[]
    for block in blocks
        T=Nemo.matrix(Nemo.ZZ,block.T);O=Nemo.matrix(Nemo.ZZ,block.O)
        X=T^(n-2)
        ub=BigInt(tr(X*O));vb=BigInt(tr(X*T*T))
        u+=ub;v+=vb
        push!(records,Dict("delta"=>block.delta,"dimension"=>size(block.T,1),
            "U"=>string(ub),"V"=>string(vb)))
        println("Exact block ",block.delta," complete in ",time()-started," seconds")
        flush(stdout)
    end
    @assert v>0 && -700741v-500000u>0
    write_result(joinpath(@__DIR__,"trial_scan.json"),Dict(
        "metadata"=>metadata("other trial lengths"),
        "floating"=>floating,"exact"=>Dict("length"=>n,"U"=>string(u),"V"=>string(v),
        "positive_margin"=>string(-700741v-500000u),"seconds"=>time()-started,
        "blocks"=>records,"scope"=>"fresh exact Julia contraction, not a Lean check")))
    println("Exact trial N=$n energy=",Float64(u//v)," seconds=",time()-started)
end
trial_main()
