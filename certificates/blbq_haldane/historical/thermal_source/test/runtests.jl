using Test, LinearAlgebra, ThermalMPOGap, ThermalMPOGap.Results, JSON3
using ThermalMPOGap.ExactDiagonalization
include("../sectors.jl")
BLAS.set_num_threads(1)
@testset "RP coefficient identities, exact rational arithmetic" begin
    r=exact_rp()
    @test -r.X==r.A
    @test -r.X^2 == -r.A/2-r.B
    # Orthogonal antisymmetric and symmetric operator spans prove the iff condition.
    @test all(dot(g,q)==0 for g in r.G,q in r.Q)
    independent=[r.G;[r.Q[a,b] for (a,b) in ((1,1),(2,2),(3,3),(1,2),(1,3),(2,3))]]
    @test det(hcat(vec.(independent)...)) != 0
end
@testset "Full dense spectrum agrees with unreduced ED" begin
    X=sz_bond()
    for L in 2:5, theta in (0.0,-0.3,-0.7)
        h=cos(theta)*X+sin(theta)*X^2
        c=Chain(3,ring_bonds(h);charge=[-1,0,1])
        direct=spectrum(c,L;dense_max=10000).E
        E=Float64[];counted=0
        for M in 0:L
            o=orbits(L,M)
            for k in 0:L÷2
                count(p->mod(k*p,L)==0,o.periods)==0 && continue
                mult=(M==0 ? 1 : 2)*((k==0 || 2k==L) ? 1 : 2)
                for b in reduced_components(L,M,k,o)
                    ev=eigvals(Symmetric(Matrix(cos(theta)*b.A+sin(theta)*b.B)))
                    for _ in 1:mult;append!(E,ev);end
                    counted+=mult*length(ev)
                end
            end
        end
        @test counted==3^L
        @test sort(E)≈direct atol=1e-10 rtol=1e-11
    end
end
qparse(s)=let a=split(s,"//");parse(BigInt,a[1])//parse(BigInt,a[2]);end
@testset "Stored rejected rational inputs" begin
    path=joinpath(@__DIR__,"../results/scaling.json")
    @test isfile(path)
    if isfile(path)
        data=JSON3.read(read(path,String))
        for row in data.best_even
            got=certified_gap_bound(qparse(row.epsilon_s_rational),
                qparse(row.epsilon_t_rational),qparse(row.beta_rational)).bound
            @test got===nothing
            @test row.bound===nothing
        end
        @test data.certified_gap===nothing
        for (name,hash) in pairs(data.input_hashes)
            @test sha256_file(joinpath(@__DIR__,"../results",string(name)))==hash
        end
        for f in filter(x->startswith(x,"spectrum_L"),readdir(joinpath(@__DIR__,"../results")))
            d=JSON3.read(read(joinpath(@__DIR__,"../results",f),String))
            @test sum(s.dimension*s.multiplicity for s in d.sectors)==3^d.L
            @test all(length(s.energies)==3^d.L for s in d.spectra)
            for spectrum in d.spectra
                c,s=cos(spectrum.theta),sin(spectrum.theta)
                mean=d.L*4s/3
                variance=d.L*(4c^2/3-4c*s/3+8s^2/9)
                d.L==2 && (variance*=2) # the two directed bonds coincide
                @test sum(spectrum.energies)/3^d.L ≈ mean atol=1e-10
                @test sum(abs2,spectrum.energies)/3^d.L ≈ variance+mean^2 rtol=1e-10
            end
        end
    end
end
