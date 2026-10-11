using Spin2Gap, Test, LinearAlgebra, Random
const SV = Spin2Gap.Variational
BLAS.set_num_threads(1)

@testset "Reduced energy and analytic gradients" begin
    for spins in ([1], [1, 2], [0, 1, 1, 2, 2, 3]), seed in (42, 71)
        model = SV.model(spins)
        rng = MersenneTwister(seed)
        x = randn(rng, length(model.t.pairs))
        x /= norm(x)
        e, g = SV.energy_gradient(model, x)
        ef, gf = SV.energy_gradient(model, x; symmetric = false)
        @test e ≈ ef atol=2e-12
        @test g ≈ gf atol=2e-10
        embedding = hcat(
            [
                SV.symmetric_expansion(v, model.symmetric_columns, length(model.t.singlets)) for
                v in eachcol(
                    Matrix{Float64}(
                        I,
                        length(model.symmetric_columns),
                        length(model.symmetric_columns),
                    ),
                )
            ]...,
        )
        @test transpose(embedding)*embedding ≈ I atol=1e-14
        reference = SV.polynomial_model(spins)
        wr, br = SV.channels(reference, x)
        wf, bf = SV.channels(model, x)
        @test wf*embedding ≈ embedding*SV.symmetric_restriction(wf, model.symmetric_columns) atol=1e-13
        @test bf*embedding ≈ embedding*SV.symmetric_restriction(bf, model.symmetric_columns) atol=1e-13
        @test wf ≈ wr atol=2e-13
        @test bf ≈ br atol=2e-13
        ep, gp = SV.energy_gradient(reference, x)
        @test e ≈ ep atol=2e-12
        @test g ≈ gp atol=2e-10
        @test e ≈ SV.ST.energy(model.t, x) atol=2e-12
        @test abs(dot(x, g)) < 1e-9
        e2, g2 = SV.energy_gradient(model, 2x)
        @test e2 ≈ e atol=2e-12
        @test 2g2 ≈ g atol=1e-9
        for _ in 1:3
            direction = randn(rng, length(x))
            direction /= norm(direction)
            h = 1e-4
            f(a) = SV.ST.energy(model.t, x + a * direction)
            difference = (f(-2h) - 8f(-h) + 8f(h) - f(2h)) / (12h)
            @test dot(g, direction) ≈ difference atol=2e-8 rtol=2e-8
        end
    end
end

@testset "Optimization and embedding" begin
    small = SV.model([1, 2])
    x = randn(MersenneTwister(81), length(small.t.pairs))
    initial, _ = SV.energy_gradient(small, x)
    best = SV.optimize(small, x; iterations = 300)
    @test best.energy < initial - 0.1
    @test best.energy ≈ SV.ST.energy(small.t, best.x) atol=2e-12
    @test best.energy < -4.74
    limited = SV.optimize_lbfgs(small, x; iterations = 500)
    @test limited.energy < -4.74
    @test limited.energy ≈ SV.ST.energy(small.t, limited.x) atol=2e-12
    metric = SV.transfer_preconditioner(small, x)
    @test dot(x, metric(x)) > 0
    preconditioned = SV.optimize_lbfgs(small, best.x+0.005x; iterations = 500, precondition = true)
    @test preconditioned.energy < -4.74
    @test preconditioned.energy ≈ SV.ST.energy(small.t, preconditioned.x) atol=2e-12
    larger = SV.model([0, 1, 1, 2, 2, 3])
    embedded = SV.embed(small.t, best.x, larger.t; noise = 0.0)
    @test first(SV.energy_gradient(larger, embedded)) ≈ best.energy atol=2e-11
    @test SV.ST.energy(larger.t, embedded) ≈ best.energy atol=2e-11
    # Exactly unused multiplets must not make a numerical search metric singular.
    metric = SV.transfer_preconditioner(larger, embedded)
    direction = randn(MersenneTwister(117), length(embedded))
    @test all(isfinite, metric(direction))
    @test dot(direction, metric(direction)) > 0
end
