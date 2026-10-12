using Test
include("../common.jl")
@testset "exact RP and rotation" begin
    for h in FIELDS; @test exact_rp(h);end
end
@testset "independent small-ring spectra and sector counts" begin
    for L in (2,4,6)
        a=all_spectra(L,3//4;verified=true)
        b=spectrum(ThermalMPOGap.Models.staggered_heisenberg(3//4).chain,L;dense_max=1000)
        @test a.spec.N==2^L
        @test maximum(abs.(a.spec.E-b.E))<1e-12
        @test sum(r["dimension"]*r["multiplicity"] for r in a.records)==2^L
        for beta in (1//2,3//1)
            z=enclosed_logz(a.records,beta)
            @test inf(z)<=log_partition(b.E,beta)<=sup(z)
        end
    end
end
@testset "stored rational seeds and bound" begin
    files=filter(f->startswith(f,"certificate_") && endswith(f,".json"),readdir(joinpath(ROOT,"results")))
    @test !isempty(files)
    for file in files
        c=JSON3.read(read(joinpath(ROOT,"results",file),String),Dict{String,Any})
        @test exact_rp(rat(c["h"]))
        for (key,L) in (("small_sectors",c["n0"]),("large_sectors",2c["n0"]))
            records=c[key]
            @test all(r["dense"] && length(r["eigenvalues"])==r["dimension"] for r in records)
            @test sum(r["dimension"]*r["multiplicity"] for r in records)==2^L
            @test all(0<=r["orthogonality_upper"]<1 && r["radius"]>0 for r in records)
        end
        beta=rat(c["beta"]);s=rat(c["epsilon_s"]);t=rat(c["epsilon_t"])
        seeds=enclosed_seeds(c["small_sectors"],c["large_sectors"],beta)
        @test seeds.s==s && seeds.t==t
        bound=certified_gap_bound(s,t,beta)
        @test bound.bound==c["bound"]
        @test bound.bound>0
    end
end
