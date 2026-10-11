using Spin2Gap, Test
import Nemo
const TB = Spin2Gap.ThermalBounds

# Independently diagonalize the four-site ring analytically:
# H=(S1+S3).(S2+S4), with sublattice spins a,b in 0:4.
function four_site_partition(beta, RR)
    value = RR(0)
    for a in 0:4, b in 0:4, j in abs(a-b):(a+b)
        e = (j*(j+1)-a*(a+1)-b*(b+1))//2
        value += (2j+1)*exp(-RR(beta)*RR(e))
    end
    return value
end

@testset "Exact rational physical projectors" begin
    bond = TB.rational_bond()
    @test [Nemo.tr(p) for p in bond.projectors] == [1, 3, 5, 7, 9]
    @test Nemo.tr(bond.h) == 0
end

@testset "One-sided Trotter partition certificates" begin
    RR = Nemo.ArbField(192)
    @test four_site_partition(0, RR) == 625
    for beta in (1//2, 2//1)
        exact = four_site_partition(beta, RR)
        previous = nothing
        for steps in (1, 2, 4, 8, 16, 64)
            bound = TB.trotter_upper(4, beta, steps)
            @test bound.trotter_trace > exact
            @test sum((r.magnetization == 0 ? 1 : 2)*r.dimension for r in bound.sectors) == 625
            isnothing(previous) || @test previous > bound.trotter_trace
            previous = bound.trotter_trace
        end
        @test previous/exact < 101//100
    end
    shifted = TB.trotter_upper(4, 1//2, 8; shift = 19//4)
    plain = TB.trotter_upper(4, 1//2, 8)
    @test Nemo.overlaps(shifted.trotter_trace, plain.trotter_trace*exp(-RR(19//2)))
    @test_throws ArgumentError TB.trotter_upper(5, 1//2, 8)
    @test_throws ArgumentError TB.trotter_upper(4, 0.5, 8)
end
