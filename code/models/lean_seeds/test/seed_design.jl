# Study `seed_design`: exact replay, forged chains and filters, cached-spectrum controls.
module TestSeedDesign
using Test, JSON3
using Main.LeanSeeds.SeedDesign
using CertifiedGaps.ExactArithmetic
@testset "exact scalar replay and forged witnesses" begin
    @test check_seed_design()
    cert=JSON3.read(read(joinpath(OUT,"certificate.json"),String),Dict{String,Any})
    for row in cert["designs"]
        input=row["input"];c=row["chain"]
        @test chain_valid(c)
        forged=deepcopy(c);forged["r"]="10//1"
        @test !chain_valid(forged)
        forged=deepcopy(c);forged["s"][2]="0//1"
        @test !chain_valid(forged)
        forged=deepcopy(c);forged["epsilon"]="1//1"
        @test !chain_valid(forged)
        for f in input["filters"]
            isnothing(f) && continue
            P=parse_rational.(f["P"]);B=parse_rational(f["B"])
            cap=parse_rational(f["cap"])
            @test filter_valid(P,B,cap)
            @test !filter_valid(P,2B,cap)
            @test !filter_valid(P,B,cap/2)
            @test all(iszero(P[i]) for i in 2:2:length(P)) # Only even moments needed.
        end
        eps=parse_rational(input["eps"])
        moments,errors=exact_sector_data(input["source"],input["mode"],eps)
        @test errors[end]==eps/2
        @test all(parse_rational(m)>0 for m in input["masses"])
    end
    @test_throws ErrorException parse_rational(0.1)
    @test capped_upper(frac(3,2),frac(1),frac(2))==frac(5,4)
    lo,hi=power_interval(frac(2),frac(5,2))
    @test lo^2<=32<=hi^2
end
@testset "Float64 controls remain explicitly non-certifying" begin
    for n in 4:10,g in 0:2
        path=joinpath(CACHE,"spectrum-$n-$g.json")
        isfile(path) || continue
        s=JSON3.read(read(path,String),Dict{String,Any})
        @test s["certified"]===false
        @test sum(x["dim"]*x["multiplicity"] for x in s["sectors"])==3^n
        @test all(issorted(x["eigenvalues"]) for x in s["sectors"])
        if n==4 && g==0
            @test abs(minimum(first(x["eigenvalues"]) for x in s["sectors"])+6)<1e-12
        end
    end
end
end
