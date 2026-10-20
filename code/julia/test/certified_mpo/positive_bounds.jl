include(joinpath(@__DIR__, "setup.jl"))
using Test, Random, LinearAlgebra
import Nemo
const PB = CM.PositiveBounds
const F = PB.F

@testset "Exact positive compression envelope" begin
    rng = MersenneTwister(55)
    for shape in ((7, 9), (10, 6))
        M = F.(rand(rng, 0:20, shape...)) ./ 20
        for rank in (2, min(shape...)), underweight in (1.0, 256.0)
            c = PB.upper_factorization(
                M,
                rank;
                iterations = 100,
                underweight,
                row_metric = collect(1:shape[1]),
                col_metric = collect(1:shape[2]),
            )
            @test all(c.U .>= 0) && all(c.V .>= 0)
            @test all(c.bound .>= M)
            @test c.rank <= rank+1
            l, r = big.(rand(rng, 0:10, shape[1])), big.(rand(rng, 0:10, shape[2]))
            @test dot(l, c.bound*r) >= dot(l, M*r)
            @test sum(abs2, c.bound) >= sum(abs2, M)
            rank == min(shape...) && @test c.bound == M
        end
    end
    z = PB.upper_factorization(zeros(F, 5, 5), 2)
    @test iszero(z.bound)
    @test_throws ArgumentError PB.upper_factorization(F[1 -1; 2 3], 1)
    RR = Nemo.ArbField(192)
    for x in (RR(0), RR(1)/3, exp(RR(-30)), exp(RR(200)))
        lo, hi = PB.rational_bracket(x)
        @test x > lo && x < hi
    end
end

@testset "Numerical covering candidates" begin
    V = [1.0 0.0 1.0; 0.0 1.0 1.0]
    b = [1.0, 2.0, 4.0]
    u = PB.covering_candidate(V, b)
    @test u ≈ [2.0, 2.0] atol=1e-10
    V = [1.0 2.0 0.0; 0.0 1.0 3.0]
    b = [0.5, 3.0, 3.0]
    @test PB.covering_candidate(V, b) ≈ [1.0, 1.0] atol=1e-10
    rng = MersenneTwister(97)
    for _ in 1:8
        V = rand(rng, 2, 4)
        b = rand(rng, 4)
        cost = V*b
        constraints = vcat(transpose(V), Matrix{Float64}(I, 2, 2))
        rhs = vcat(b, zeros(2))
        candidates = Vector{Float64}[]
        for i in 1:5, j in (i+1):6
            A = constraints[[i, j], :]
            abs(det(A)) > 1e-12 || continue
            candidate = A\rhs[[i, j]]
            all(constraints*candidate .>= rhs .- 1e-10) && push!(candidates, candidate)
        end
        optimum = minimum(dot(cost, x) for x in candidates)
        u = PB.covering_candidate(V, b)
        @test all(transpose(V)*u .>= b .- 1e-9)
        @test dot(cost, u) ≈ optimum atol=1e-9
        scales = [1e-8, 1e8]
        rescaled = PB.covering_candidate(V .* scales, b) .* scales
        @test all(transpose(V)*rescaled .>= b .- 1e-9)
        @test dot(cost, rescaled) ≈ optimum atol=1e-9
    end
    M = F.(rand(rng, 1:20, 8, 9)) ./ 20
    c = PB.upper_factorization(M, 3; iterations = 100, cover_sweeps = 2)
    @test all(c.bound .>= M)
    @test all(c.U .>= 0) && all(c.V .>= 0)
    split = PB.upper_factorization(M, 3; iterations = 100, correction_rank = 4)
    @test all(split.bound .>= M)
    @test all(split.U .>= 0) && all(split.V .>= 0)
    @test split.rank <= 7
    baseline = PB.upper_factorization(M, 3; iterations = 100)
    @test sum(abs2, split.bound) <= sum(abs2, baseline.bound)
    @test sum(abs2, c.bound) <= sum(abs2, baseline.bound)
end
