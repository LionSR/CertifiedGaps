# Study `spectrum`: the Int128 floor primitive against BigInt, fixed-point against Float64
# filtering, then the data consistency checks (about 15 s in total).
module TestSpectrum
using Test, LinearAlgebra, Random
using Main.LeanSeeds.Spectrum

@testset "Integer recurrence primitives" begin
    rng=MersenneTwister(72)
    for _ in 1:1000
        x=(Int128(rand(rng,Int64))<<48)+rand(rng,UInt32)
        c=Int128(rand(rng,-(Int64(1)<<43):(Int64(1)<<43)))
        @test BigInt(mulfloor(x,c))==fld(BigInt(x)*BigInt(c),BigInt(1)<<40)
    end
    for twist in TWISTS
        H,_=hamiltonian(6,0,twist)
        A,B=integer_hamiltonian(H)
        @test isnothing(B) ? A≈H : A+cis(2pi/3)*B≈H
        Q,_=product_seeds(H,4,123)
        coeffs=cheb_coefficients(-9.,-7.,6.,16)
        Y,_=float_filter(H,Q,coeffs)
        for bits in (64,128)
            Z,_=fixed_filter((A,B),Q,coeffs,bits)
            @test norm(Y-Z)<1e-9
        end
    end
end
check_hamiltonian()
check_data()
end
