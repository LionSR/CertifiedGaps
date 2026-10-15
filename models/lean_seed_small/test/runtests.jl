using Test, JSON
include("../check.jl")

@testset "exact certificate and mutation rejection" begin
    c=JSON.parsefile(joinpath(OUT,"certificate.json"))
    @test verify_certificate(c)
    bad=deepcopy(c); bad["points"][1]["lambda"][1]="0//1"
    @test_throws AssertionError verify_certificate(bad)
    bad=deepcopy(c)
    bad["points"][2]["anchor_20H"]=exactstr(exactint(bad["points"][2]["anchor_20H"])+1)
    @test_throws AssertionError verify_certificate(bad)
    bad=deepcopy(c); bad["points"][1]["packed"][1]="0//1"
    @test_throws AssertionError verify_certificate(bad)
    bad=deepcopy(c); bad["points"][1]["residual_sq"][1]="0//1"
    @test_throws AssertionError verify_certificate(bad)
end

@testset "Hamiltonian boundary convention" begin
    # For length two both directed periodic bonds occur, as in RPRing.hamiltonian.
    @test diagonal(2,0)==-42
    @test diagonal(2,1)==40
    @test apply_h(2,BigInt[1,0,0,0])==BigInt[-42,-5,-5,0]
end
