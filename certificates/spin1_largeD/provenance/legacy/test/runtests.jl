using Test
include("../seed.jl")

@testset "Stored spin-one large-D certificates" begin
    @test VERSION==v"1.12.4"
    @test BLAS.get_num_threads()==1
    @test sha256_file(joinpath(root,"PROMPT.txt")) ==
        "00d49b8627a72d4ce4841d0cfe85bb0accc9f5a9d13febb9b124a2c488cbb975"
    names=filter(x->startswith(x,"certificate_D")&&endswith(x,".json"),
        readdir(joinpath(root,"results")))
    @test length(names)>=2
    for name in names
        data=JSON3.read(read(joinpath(root,"results",name),String),Dict{String,Any})
        D=parse(Rational{BigInt},data["D"])
        beta=parse(Rational{BigInt},data["beta"])
        es=parse(Rational{BigInt},data["eps_s"])
        et=parse(Rational{BigInt},data["eps_t"])
        n=Int(data["n0"])
        @test rp_identity(D)
        @test iseven(n) && data["certified"]
        @test 0<=es<1 && 0<=et<1 && beta>0
        for (file,hash) in data["source_hashes"]
            @test sha256_file(joinpath(root,file))==hash
        end
        for (file,hash) in data["core_hashes"]
            @test sha256_file(joinpath(root,"../../julia",file))==hash
        end
        for (file,hash) in data["enclosure_hashes"]
            @test sha256_file(joinpath(root,file))==hash
        end
        a=JSON3.read(read(joinpath(root,"results","enclosed_D$(numerator(D))_$(denominator(D))_L$n.json"),String),Dict{String,Any})
        b=JSON3.read(read(joinpath(root,"results","enclosed_D$(numerator(D))_$(denominator(D))_L$(2n).json"),String),Dict{String,Any})
        for spec in (a,b)
            @test spec["complete"] && spec["rigorous"]
            @test sum(s["dimension"]*s["multiplicity"] for s in spec["sectors"])==3^spec["L"]
            @test length(spec["energies"])==length(spec["radii"])==3^spec["L"]
            @test all(s["levels"]==s["dimension"] for s in spec["sectors"])
            @test all(s["orthogonality_upper"]<1 for s in spec["sectors"])
            for (file,hash) in spec["source_hashes"]
                @test sha256_file(joinpath(root,file))==hash
            end
        end
        za=enclosed_logz(a,beta); zb=enclosed_logz(b,beta); zbb=enclosed_logz(b,2beta)
        @test es>=Rational{BigInt}(sup(1-exp(zb-2za)))
        @test et>=Rational{BigInt}(sup(1-exp(zbb-2zb)))
        result=certified_gap_bound(es,et,beta)
        @test result.bound!==nothing
        @test result.bound==data["bound"]>0
        println("D=$D beta=$beta uniform bound=$(result.bound), even L >= $n sites")
    end
end
