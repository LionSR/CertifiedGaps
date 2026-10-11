include(joinpath(@__DIR__, "setup.jl"))
using Test
import Nemo
const PG = CM.PolynomialBounds

@testset "Polynomial thermal majorant" begin
    RR = Nemo.ArbField(192)
    for power in (4, 12, 28)
        c = PG.seed(4, 4, power)
        @test c.delta*2c.m == c.beta
        @test RR(c.remainder) > exp(RR(c.delta*c.spectral_radius))*RR(c.delta*c.spectral_radius)^3/6
        for k in -20:20
            q = RR(PG.polynomial(c, k))
            exact = exp(-RR(c.delta)*k)
            @test q > exact
            @test q-exact < 2RR(c.remainder)
        end
        # The four-site ring is (S1+S3).(S2+S4), so its full spectrum follows
        # from adding two pairs of spin-2 representations independently.
        spectrum = [
            (div(j*(j+1)-a*(a+1)-b*(b+1), 2), 2j+1) for a in 0:4 for b in 0:4 for
            j in abs(a-b):(a+b)
        ]
        @test sum(last, spectrum) == 625
        exact_trace = sum(degeneracy*exp(-RR(c.beta)*e) for (e, degeneracy) in spectrum)
        polynomial_trace =
            exp(RR(c.beta)*4)*sum(
                degeneracy*RR(PG.polynomial(c, e+4))^(2c.m) for (e, degeneracy) in spectrum
            )
        @test polynomial_trace > exact_trace
        @test polynomial_trace < RR(c.inflation_upper)*exact_trace
    end
    target = PG.seed(256, 64, 28)
    @test target.inflation_upper < PG.F(1001, 1000)
    @test PG.seed(256, 64, 30).inflation_upper < target.inflation_upper
end

@testset "Exact positive Taylor sum of squares" begin
    quartic = PG.even_taylor_sos(4)
    @test quartic.weights == [1, 1//12, 1//72]
    @test quartic.polynomials == [[1, 1//2, 1//12], [0, 1, 1//2], [0, 0, 1]]
    for degree in (2, 4, 6, 8, 12, 16, 24, 32)
        sos = PG.even_taylor_sos(degree)
        @test all(>(0), sos.weights)
        for x in (-7//3, -1//10, 0//1, 3//5, 5//2)
            expected = sum(PG.F(x)^j/factorial(big(j)) for j in 0:degree)
            actual = sum(
                w*sum(c*PG.F(x)^(j-1) for (j, c) in enumerate(p))^2 for
                (w, p) in zip(sos.weights, sos.polynomials)
            )
            @test actual == expected
            @test actual > 0
        end
    end
end

@testset "Higher-order majorants against the complete four-site spectrum" begin
    RR = Nemo.ArbField(384)
    spectrum =
        [(div(j*(j+1)-a*(a+1)-b*(b+1), 2), 2j+1) for a in 0:4 for b in 0:4 for j in abs(a-b):(a+b)]
    for degree in (1, 4, 8, 12), power in (4, 8)
        c = PG.taylor_seed(4, 4, power, degree)
        for k in -20:20
            value = PG.taylor_polynomial(c, k)
            exact = exp(-RR(c.delta)*k)
            @test RR(value) > exact
            @test RR(value)-exact < 2RR(c.remainder)
            @test abs(value) <= c.operator_upper
        end
        exact_trace = sum(degeneracy*exp(-RR(c.beta)*e) for (e, degeneracy) in spectrum)
        polynomial_trace =
            exp(RR(c.beta)*4)*sum(
                degeneracy*RR(PG.taylor_polynomial(c, e+4))^(2c.m) for (e, degeneracy) in spectrum
            )
        @test polynomial_trace > exact_trace
        @test polynomial_trace < exp(RR(c.log_inflation_upper))*exact_trace
    end
end
