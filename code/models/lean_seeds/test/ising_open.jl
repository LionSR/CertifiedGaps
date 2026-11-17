# Study `ising_open`: exact certificate replay, mutation rejection and the free-end convention.
module TestIsingOpen
using Test, JSON
using Main.LeanSeeds.IsingOpenSeed
using Main.LeanSeeds.IsingSeed: exactstr, exactint

@testset "exact certificate and mutation rejection" begin
    c = JSON.parsefile(joinpath(OUT, "certificate.json"))
    @test verify_certificate(c)
    bad = deepcopy(c); bad["points"][1]["lambda"][1] = "0//1"
    @test_throws AssertionError verify_certificate(bad)
    bad = deepcopy(c); bad["points"][2]["packed"][1] = "0//1"
    @test_throws AssertionError verify_certificate(bad)
    # rho is one below the largest value the weights allow, so rho + 2 must fail
    bad = deepcopy(c)
    bad["points"][3]["rho_20H"] = exactstr(exactint(bad["points"][3]["rho_20H"]) + 2)
    @test_throws AssertionError verify_certificate(bad)
    bad = deepcopy(c); bad["points"][3]["weights"] = "0//1"
    @test_throws AssertionError verify_certificate(bad)
end

@testset "free-end convention" begin
    # The chain of two sites has one bond, unlike the ring of two sites.
    @test diagonal(2, 0) == -22
    @test diagonal(2, 1) == 20
    @test apply_h(2, BigInt[1, 0, 0, 0]) == BigInt[-22, -5, -5, 0]
    # flipping an end spin breaks one bond (cost 40 in 20 H), a bulk spin two
    @test diagonal(8, 1) - diagonal(8, 0) == 42
    @test diagonal(8, 2) - diagonal(8, 0) == 82
end
end
