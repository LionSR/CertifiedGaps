include("ladder.jl")
using Test
@test exact_rp(4//1)
for L in (2,3,4,5)
    J=7//2
    # Independent bit-basis full matrix, without any symmetry reduction.
    H=zeros(4^L,4^L)
    for s in 0:4^L-1, (t,h) in transitions(s,L,J); H[t+1,s+1]=Float64(h); end
    x=[0.0 0.5;0.5 0.0]; a=[0.0 0.5;-0.5 0.0]; z=[0.5 0;0 -0.5]
    h=kron(x,x)-kron(a,a)+kron(z,z)
    bonds=[(2j+1,2j+2,Float64(J)*h) for j in 0:L-1]
    append!(bonds,[(2j+leg,2mod(j+1,L)+leg,h) for j in 0:L-1 for leg in 1:2])
    chain=Chain(2,_ -> bonds)
    @test H == Matrix(sector_hamiltonian(chain,2L,collect(0:4^L-1)))
    full,fullradius,_=verified_eigen(H,0.0)
    weight2=L*Float64(J)^2 + (L==2 ? 8 : 2L)
    @test sum(abs2,full)/4^L ≈ 3weight2/16 atol=1e-10
    red=full_spectrum(L,J)
    @test red.energies ≈ full atol=1e-11
    checked=full_spectrum(L,J;verified=true)
    @test all(abs.(checked.energies-full) .<= checked.radii .+ fullradius)
end
println("Symmetry spectra and eigenvalue enclosures match independent full matrices.")
