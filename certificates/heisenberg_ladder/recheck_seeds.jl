include("certify.jl")
using Test
function recheck()
    c=JSON3.read(read(joinpath(@__DIR__,"results/certificate.json"),String))
    parseq(s)=parse(Rational{BigInt},s)
    specs=Dict()
    for L in (Int(c.n0),2Int(c.n0))
        s=JSON3.read(read(joinpath(@__DIR__,"results/verified_L$(L).json"),String))
        specs[L]=(;energies=Float64.(parseq.(s.energies)),radii=Float64.(parseq.(s.radii)))
    end
    beta=parseq(c.beta)
    a=enclosed_logZ(specs[Int(c.n0)],beta)
    b=enclosed_logZ(specs[2Int(c.n0)],beta)
    d=enclosed_logZ(specs[2Int(c.n0)],2beta)
    @test Rational{BigInt}(sup(-expm1(b-2a))) <= parseq(c.eps_s)
    @test Rational{BigInt}(sup(-expm1(d-2b))) <= parseq(c.eps_t)
    println("Stored complete eigenvalue enclosures reproduce both rational seed upper bounds.")
end
recheck()
