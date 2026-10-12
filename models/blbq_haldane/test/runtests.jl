using Test, LinearAlgebra, CertifiedGaps, CertifiedGaps.Results, CertifiedGaps.Seeds, JSON3
using CertifiedGaps.ExactDiagonalization
include("../sectors.jl")
include("../result_io.jl")
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

include("su2_small.jl")
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
const RESULTS=joinpath(@__DIR__,"../results")
const THETAS=[0.0,-0.1,-0.3,-0.5,-0.6,-0.7]
@testset "Stored rejected rational inputs and summary tables" begin
    path=joinpath(RESULTS,"scaling.json")
    @test isfile(path)
    data=JSON3.read(read(path,String))
    for row in data.best_even
        got=certified_gap_bound(qparse(row.epsilon_s_rational),
            qparse(row.epsilon_t_rational),qparse(row.beta_rational)).bound
        @test got===nothing
        @test row.bound===nothing
    end
    @test data.certified_gap===nothing
    @test all(r->r.bound===nothing,data.grid_maxima)
    # Complete spectra (spectrum_L*.json[.gz]) were removed from the repository to keep the
    # negative result small; their hashes stay in scaling.json and archive.json. Files that
    # are still present must match their recorded hashes.
    for (name,hash) in pairs(data.input_hashes)
        file=joinpath(RESULTS,string(name))
        if startswith(string(name),"spectrum_L")
            @test !isfile(file) || sha256_file(file)==hash
        else
            @test sha256_file(file)==hash
        end
    end
    # The L14 multiplet-resolved ground energies and gaps must agree with the independent
    # low-energy reference (numerical consistency only, no enclosure claim).
    reference=JSON3.read(read(joinpath(RESULTS,"reference_L14.json"),String))
    for (t,theta) in enumerate(THETAS)
        row=only(filter(r->r.theta==theta && r.n==7,data.energy_scaling))
        @test row.complete_spectrum_2n
        @test row.ground_2n ≈ reference.summary[t].ground atol=1e-8
        @test row.gap_2n ≈ reference.summary[t].gap atol=1e-8
    end
end

@testset "Small-ring rows of the summary tables reproduced from scratch" begin
    data=JSON3.read(read(joinpath(RESULTS,"scaling.json"),String))
    X=sz_bond()
    for (t,theta) in enumerate(THETAS)
        h=cos(theta)*X+sin(theta)*X^2
        specs=Dict(L=>spectrum(Chain(3,ring_bonds(h);charge=[-1,0,1]),L;dense_max=1000)
            for L in 2:6)
        for (L,s) in specs
            # Analytic infinite-temperature moments of the complete spectrum.
            c,si=cos(theta),sin(theta)
            mean=L*4si/3
            variance=L*(4c^2/3-4c*si/3+8si^2/9)
            L==2 && (variance*=2) # the two directed bonds coincide
            @test sum(s.E)/3^L ≈ mean atol=1e-10
            @test sum(abs2,s.E)/3^L ≈ variance+mean^2 rtol=1e-10
        end
        results=scan_seeds(spectral_defects(specs),2:3,0.125:0.125:64.0)
        for n in 2:3
            best=argmax(r->min(1-r.εₛ,1-r.εₜ),filter(r->r.n₀==n,results))
            row=only(filter(r->r.theta==theta && r.n==n,data.grid_maxima))
            @test best.β₀==row.beta
            @test 1-best.εₛ ≈ row.spatial rtol=1e-9
            @test 1-best.εₜ ≈ row.thermal rtol=1e-9
            @test best.bound===nothing
            srow=only(filter(r->r.theta==theta && r.n==n,data.energy_scaling))
            @test specs[2n].E[1]-2specs[n].E[1] ≈ srow.delta atol=1e-9
        end
    end
end

@testset "Complete weighted-spectrum partition function" begin
    weighted=WeightedSpectrum([0.0,1.0,2.0],[1,3,5],9,true)
    expanded=PartialSpectrum([0.0;ones(3);fill(2.0,5)],9,true)
    for beta in (0.125,1.0,8.0)
        @test logZ_bounds(weighted,beta)[1]≈logZ_bounds(expanded,beta)[1] atol=1e-14
    end
    @test spectrum_moment(weighted,1)==spectrum_moment(expanded,1)
    @test spectrum_moment(weighted,2)==spectrum_moment(expanded,2)
end
