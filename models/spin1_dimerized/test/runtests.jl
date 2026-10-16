# Fast checker: exact RP and cell reflection, the symmetry decomposition against direct ED,
# exact larger-period projectors and the verified block path on a small ring, the stored seeds
# and bounds with their log-partition enclosures and provenance, and a fresh recertification
# of the delta = 9/10 seed (rings of 4 and 8 spins) checked against the stored enclosures.
using Test
using SHA
include(joinpath(@__DIR__, "..", "src", "common.jl"))
include(joinpath(@__DIR__, "..", "src", "rp.jl"))
include(joinpath(@__DIR__, "..", "src", "enclosure.jl"))
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
        @test check_rp(parse_rational(data["delta"]))
        n=parse(Int,data["n0"]);@test n in (2,3)
        @test data["positive_transfer_gram_proof"]
        es=parse_rational(data["epsilon_s"]);et=parse_rational(data["epsilon_t"])
        beta=parse_rational(data["beta0"])
        bound=certified_gap_bound(es,et,beta).bound
        @test bound!==nothing && bound>0
        @test bound==data["bound"]
        @test Rational{BigInt}(bound)==parse_rational(data["bound_exact_binary"])
        ti=CertifiedGaps.Recurrence.to_interval
        z1hi=ti(parse_rational(data["logZ_n_beta"]["upper"]))
        z2lo=ti(parse_rational(data["logZ_2n_beta"]["lower"]))
        z2hi=ti(parse_rational(data["logZ_2n_beta"]["upper"]))
        z22lo=ti(parse_rational(data["logZ_2n_2beta"]["lower"]))
        @test Rational{BigInt}(sup(1-exp(z2lo-2z1hi)))<=es
        @test Rational{BigInt}(sup(1-exp(z22lo-2z2hi)))<=et
        for (path,hash) in data["metadata"]["source_sha256"]
            if haskey(data["metadata"], "ThermalMPOGap")
                # The stored hashes are provenance of the original run, made before the
                # core was renamed and the directory moved from certificates/ to models/.
                # They are checked against that archived commit when it is reachable
                # (full clone); a shallow clone skips this provenance check only.
                revision="6c375b5981f1f7c7813bce53ca01a1ce1c04c9f2"
                relative=relpath(normpath(joinpath(ROOT,path)),REPO)
                relative=replace(relative,r"^models/spin1_dimerized/"=>"certificates/spin1_dimerized/")
                spec="$revision:$relative"
                if success(pipeline(`git -C $REPO cat-file -e $spec`;stderr=devnull))
                    @test bytes2hex(SHA.sha256(read(`git -C $REPO show $spec`)))==hash
                else
                    @test_skip "archived source $spec unavailable"
                end
            elseif startswith(normpath(path),"..")
                # Outside this directory (the shared core): provenance only, not compared.
                continue
            else
                @test sha256_file(joinpath(ROOT,path))==hash
            end
        end
        @test sum(parse(Int,r["dimension"])*parse(Int,r["multiplicity"]) for r in data["sectors_n"])==9^n
        @test sum(parse(Int,r["dimension"])*parse(Int,r["multiplicity"]) for r in data["sectors_2n"])==9^(2n)
    end
end
@testset "Exact larger-period projectors" begin
    for cells in (3,4,6)
        states=charge_groups(cells)[end-2]
        for k in 0:fld(cells,2),reflection in (k==0 || 2k==cells ? (-1,1) : (1,))
            _,distance=verified_basis(states,cells,k,0,reflection)
            @test distance<1e-8
        end
    end
end
@testset "Verified block basis and Hamiltonian" begin
    delta=4//5;cells=2
    vals,radii,records,_=verified_block_spectrum(cells,delta)
    reference,_,_=dense_spectrum(cells,delta)
    perm=sortperm(vals)
    @test all(abs(vals[perm[i]]-reference.E[i])<Float64(radii[perm[i]]) for i in eachindex(vals))
    @test sum(r["dimension"]*r["multiplicity"] for r in records)==9^cells
end
@testset "Fresh recertification of delta = 9/10" begin
    # Rebuilds every block, eigensystem and radius of the stored delta = 9/10 certificate from
    # the model. LAPACK output differs slightly across machines, so the fresh seeds need not be
    # bit-identical; the fresh certificate must be valid on its own, and its log Z enclosures
    # must overlap the stored ones because both enclose the same exact numbers.
    stored=read_exact_json(joinpath(ROOT,"results","certificate_d9.json"))
    path=tempname()*".json"
    fresh=redirect_stdout(devnull) do
        certify_selected(9//10; n0=2, beta0=7//2, output_path=path)
    end
    rm(path; force=true)
    @test fresh["bound"]>0
    @test certified_gap_bound(parse_rational(fresh["epsilon_s"]),parse_rational(fresh["epsilon_t"]),
        7//2).bound==fresh["bound"]
    for key in ("logZ_n_beta","logZ_2n_beta","logZ_2n_2beta")
        a=parse_rational.((stored[key]["lower"],stored[key]["upper"]))
        b=parse_rational.((fresh[key]["lower"],fresh[key]["upper"]))
        @test max(a[1],b[1])<=min(a[2],b[2])
    end
    # Machine-level differences move the bound in the ninth decimal, not more.
    @test abs(fresh["bound"]-stored["bound"])<1e-8
end
