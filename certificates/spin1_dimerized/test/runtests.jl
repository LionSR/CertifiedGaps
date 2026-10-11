using Test
include("../common.jl")
include("../rp.jl")
@testset "Exact RP and cell reflection" begin
    for delta in DELTAS;@test check_rp(delta);end
end
@testset "Complete symmetry decomposition" begin
    delta=1//2;cells=2;groups=charge_groups(cells)
    reference=Float64[]
    for m in 0:2cells
        A,den=integer_hamiltonian(cells,delta,groups[m+1])
        E=eigvals(Symmetric(Matrix(A)/den))
        append!(reference,E);m>0 && append!(reference,E)
    end
    spectrum,records,_=dense_spectrum(cells,delta)
    @test length(spectrum.E)==9^cells
    @test maximum(abs,spectrum.E-sort(reference))<1e-12
    # At the dimer limit each bond has energies -2, -1, 1 with degeneracy 1,3,5.
    decoupled,_,_=dense_spectrum(cells,1//1)
    single=[-4;fill(-2,3);fill(2,5)]
    @test maximum(abs,decoupled.E-sort([a+b for a in single for b in single]))<1e-12
end
@testset "Stored rational seed and bound" begin
    files=filter(f->startswith(f,"certificate_d")&&endswith(f,".json"),readdir(joinpath(ROOT,"results")))
    @test !isempty(files)
    for file in files
        data=read_exact_json(joinpath(ROOT,"results",file))
        @test check_rp(parserat(data["delta"]))
        n=parse(Int,data["n0"]);@test n in (2,3)
        @test data["positive_transfer_gram_proof"]
        es=parserat(data["epsilon_s"]);et=parserat(data["epsilon_t"])
        beta=parserat(data["beta0"])
        bound=certified_gap_bound(es,et,beta).bound
        @test bound!==nothing && bound>0
        @test bound==data["bound"]
        @test Rational{BigInt}(bound)==parserat(data["bound_exact_binary"])
        ti=ThermalMPOGap.Recurrence.to_interval
        z1hi=ti(parserat(data["logZ_n_beta"]["upper"]))
        z2lo=ti(parserat(data["logZ_2n_beta"]["lower"]))
        z2hi=ti(parserat(data["logZ_2n_beta"]["upper"]))
        z22lo=ti(parserat(data["logZ_2n_2beta"]["lower"]))
        @test Rational{BigInt}(sup(1-exp(z2lo-2z1hi)))<=es
        @test Rational{BigInt}(sup(1-exp(z22lo-2z2hi)))<=et
        for (path,hash) in data["metadata"]["source_sha256"]
            @test sha256_file(joinpath(ROOT,path))==hash
        end
        @test sum(parse(Int,r["dimension"])*parse(Int,r["multiplicity"]) for r in data["sectors_n"])==9^n
        @test sum(parse(Int,r["dimension"])*parse(Int,r["multiplicity"]) for r in data["sectors_2n"])==9^(2n)
    end
end
@testset "Exact larger-period projectors" begin
    include("../verified_blocks.jl")
    for cells in (3,4,6)
        states=charge_groups(cells)[end-2]
        for k in 0:fld(cells,2),reflection in (k==0 || 2k==cells ? (-1,1) : (1,))
            _,distance=verified_basis(states,cells,k,0,reflection)
            @test distance<1e-8
        end
    end
end
