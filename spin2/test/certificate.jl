using Spin2Gap, Test, Random, JSON3
const SC = Spin2Gap.Certificate
const ST = SC.ST
import Nemo

@testset "Complete virtual representation bases" begin
    for j in 0:3, k in 0:3
        b = SC.representation_basis(j, k)
        @test sum(2s+1 for s in keys(b)) == (2j+1)*(2k+1)
    end
end

@testset "Reduced exact and ball contractions" begin
    for spins in ([1], [1, 2], [0, 1, 1, 2, 2, 3])
        t = ST.templates(spins)
        parameters = big.(rand(MersenneTwister(93), -10:10, length(t.pairs)))
        P = ST.polynomial_tensors(t, parameters)
        @test SC.check_tensor_covariance(t, P)
        corrupted = deepcopy(P)
        corrupted[1][1, 1] += 1
        @test_throws AssertionError SC.check_tensor_covariance(t, corrupted)
        blocks = ST.exact_blocks(t, P)
        reduced = SC.reduced_blocks(t, P; verify_images = true)
        by_covariance = SC.reduced_blocks(t, P; verify_images = false)
        @test all(a.T == b.T && a.O == b.O for (a, b) in zip(reduced, by_covariance))
        simultaneous = SC.ball_contract_many(reduced, [4, 8]; scale = 10)
        split = SC.split_transpose_rotation(reduced)
        @test sum((2b.spin+1)*length(b.pairs) for b in reduced) == t.D^2
        for n in (4, 8)
            direct = ST.exact_contract(blocks, n)
            actual = SC.exact_contract(reduced, n)
            @test actual.U == direct.U
            @test actual.V == direct.V
            split_exact = SC.exact_contract(split, n)
            @test split_exact.U == actual.U && split_exact.V == actual.V
            ball = SC.ball_contract(reduced, n; scale = 10)
            @test Nemo.contains(ball.U, Nemo.QQ(actual.U, big(10)^(2n)))
            @test Nemo.contains(ball.V, Nemo.QQ(actual.V, big(10)^(2n)))
            @test Nemo.contains(ball.energy, Nemo.QQ(actual.U, actual.V))
            @test Nemo.contains(simultaneous[n].U, Nemo.QQ(actual.U, big(10)^(2n)))
            @test Nemo.contains(simultaneous[n].V, Nemo.QQ(actual.V, big(10)^(2n)))
            @test Nemo.contains(simultaneous[n].energy, Nemo.QQ(actual.U, actual.V))
        end
    end
end

@testset "Independent saved integer witness" begin
    input = JSON3.read(
        read(joinpath(@__DIR__, "../results/trial-certificate.json"), String),
        Dict{String,Any},
    )
    t = ST.templates(Int.(input["virtual_spins"]))
    P = ST.polynomial_tensors(t, parse.(BigInt, input["integer_parameters"]))
    reduced = SC.reduced_blocks(t, P; verify_images = true)
    simultaneous = SC.ball_contract_many(reduced, [128, 256, 512]; scale = 10^8)
    @test length(reduced) == 7
    for expected in input["certificates"]
        ball = SC.ball_contract(reduced, expected["length"]; scale = 10^8)
        @test Nemo.contains(
            ball.energy,
            Nemo.QQ(parse(BigInt, expected["U"]), parse(BigInt, expected["V"])),
        )
        @test ball.energy < Nemo.QQ(-9519, 2000)
        @test Nemo.contains(
            simultaneous[expected["length"]].energy,
            Nemo.QQ(parse(BigInt, expected["U"]), parse(BigInt, expected["V"])),
        )
    end
end
