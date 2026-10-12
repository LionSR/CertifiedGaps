include("certify.jl")
using .LongitudinalCertificate, LinearAlgebra, CertifiedGaps, CertifiedGaps.Results
using CertifiedGaps.Models, CertifiedGaps.ExactDiagonalization, Test
const C=LongitudinalCertificate
start=time(); comparisons=Any[]; maxerror=0.0
@testset "Independent computational-basis and free-fermion cross-checks" begin
    Z=[1.0 0;0 -1]; X=[0.0 1;1 0]
    for L in (2,4,6), (g,h) in ((1//4,1//20),(1//1,1//4),(2//1,1//1))
        model=Chain(2,ring_bonds(-kron(Z,Z));onsite=-Float64(g)*X-Float64(h)*Z)
        direct=spectrum(model,L;dense_max=1<<L).E
        momentum=C.ring_spectrum(L,g,h).E
        @test maximum(abs.(direct-momentum))<2e-12
    end
    for L in (2,4,6,8,10,12)
        s=C.ring_spectrum(L,2//1,0//1)
        for b in (0.25,1.0,4.0,16.0,64.0)
            a=C.logz(s.E,b); f=Float64(tfim_logZ(L,b,2))
            error=abs(a-f); global maxerror=max(maxerror,error)
            @test error<2e-10
            push!(comparisons,Dict("L"=>L,"beta"=>C.str(b),"ed_logZ"=>a,
                "free_fermion_logZ"=>f,"absolute_error"=>error))
        end
    end
    # For g=0 the exact computational spectrum is rational and provides enclosure tests.
    for L in (2,4,6)
        s=C.ring_spectrum(L,0//1,1//10;rigorous=true)
        exactlevels=sort([-(L-2count_ones(xor(t,C.rotate(t,L))))-(1//10)*(L-2count_ones(t))
            for t in 0:(1<<L)-1])
        radius=maximum(C.rat(b["radius"]) for b in s.blocks)
        @test all(abs(C.exact(v)-e)<=radius for (v,e) in zip(s.E,exactlevels))
    end
end
resources=Any[]
for L in (12,16,18,20,22,24)
    # Exact dimension of k=1: all shorter-period orbits excluded. Möbius inversion.
    primitive=Dict{Int,BigInt}()
    for p in 1:L
        L%p==0 || continue
        primitive[p]=big(2)^p-sum(v for (q,v) in primitive if p%q==0;init=big(0))
    end
    d=primitive[L]÷L
    push!(resources,Dict("L"=>L,"k1_dimension"=>string(d),
        "one_complex_dense_matrix_bytes"=>string(16d*d),
        "one_real_dense_matrix_bytes"=>string(8d*d)))
end
write_result(joinpath(C.ROOT,"results","audit.json"),Dict("crosschecks"=>comparisons,
    "max_error"=>maxerror,"resources"=>resources,"seconds"=>time()-start,
    "metadata"=>C.metadata(Dict("kind"=>"numerical crosschecks and exact dimension counts"))))
