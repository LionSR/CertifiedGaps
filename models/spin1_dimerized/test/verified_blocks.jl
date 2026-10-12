using Test
include("../certify.jl")
@testset "Verified block basis and Hamiltonian" begin
    delta=4//5;cells=2
    vals,radii,records,_=verified_block_spectrum(cells,delta)
    reference,_,_=dense_spectrum(cells,delta)
    perm=sortperm(vals)
    @test all(abs(vals[perm[i]]-reference.E[i])<Float64(radii[perm[i]]) for i in eachindex(vals))
    @test sum(r["dimension"]*r["multiplicity"] for r in records)==9^cells
    # Check exact projectors also for all larger supported cell periods.
    for cells in (3,4,6)
        states=charge_groups(cells)[end-2]
        for k in 0:fld(cells,2), reflection in (k==0 || 2k==cells ? (-1,1) : (1,))
            B,d=verified_basis(states,cells,k,0,reflection)
            @test d<1e-8
        end
    end
end
