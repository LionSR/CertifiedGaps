using Test, ThermalMPOGap, ThermalMPOGap.Results
include("src/Ladder.jl")
using .Ladder
root=@__DIR__
p=joinpath(root,get(ARGS,1,"results/audit_2026-10-11_replay.json"))
c=read_exact_json(p)
started=time()
original=joinpath(root,"results","certificate_g3_1.json")
old=read_exact_json(original)
checks=@testset "Fresh spectral replay artifact" begin
    @test sha256_file(joinpath(root,"results","certificate_g3_1.json"))==c["original_certificate_sha256"]
    for (path,hash) in c["source_hashes"]
        @test sha256_file(joinpath(root,path))==hash
    end
    @test all(c["agreement"][key] for key in (
        "all_eigenvalues_identical","all_radii_identical","all_residual_bounds_identical",
        "all_orthogonality_bounds_identical","all_construction_bounds_identical",
        "epsilon_s_identical","epsilon_t_identical","bound_identical","original_bound_supported"))
    # Independently reconstruct comparison claims rather than trusting saved booleans.
    @test length(c["spectra"])==length(old["spectra"])
    for (sp,osp) in zip(c["spectra"],old["spectra"])
        @test sp["L"]==osp["L"]
        @test length(sp["sectors"])==length(osp["sectors"])
        for (sec,osec) in zip(sp["sectors"],osp["sectors"])
            @test (sec["parity"],sec["character"],sec["dimension"])==
                (osec["parity"],osec["character"],osec["dimension"])
            @test sec["eigenvalues"]==osec["eigenvalues"]
            for key in ("radius","residual","orthogonality","construction")
                @test sec["verification"][key]==osec["verification"][key]
            end
        end
    end
    @test c["epsilon_s"]==old["epsilon_s"]
    @test c["epsilon_t"]==old["epsilon_t"]
    @test c["bound"]==old["bound"]
    @test c["logZ"]==old["logZ"]
    @test c["recurrence"]==old["recurrence"]

    beta=parse(Rational{BigInt},c["beta"])
    specs=[(;sectors=[(;eigenvalues=Float64.(sec["eigenvalues"]),
        verification=(;radius=sec["verification"]["radius"])) for sec in sp["sectors"]])
        for sp in c["spectra"]]
    ln=logz_interval(specs[1],beta); l2=logz_interval(specs[2],beta)
    l22=logz_interval(specs[2],2beta)
    ls=Ladder.down(()->l2[1]-2ln[2]); lt=Ladder.down(()->l22[1]-2l2[2])
    es=Rational{BigInt}(-Ladder.down(()->expm1(ls)))
    et=Rational{BigInt}(-Ladder.down(()->expm1(lt)))
    @test es==parse(Rational{BigInt},c["epsilon_s"])
    @test et==parse(Rational{BigInt},c["epsilon_t"])
    @test certified_gap_bound(es,et,beta).bound==c["bound"]
    @test Rational{BigInt}(c["bound"])==parse(Rational{BigInt},c["bound_exact_binary"])
    for sp in c["spectra"]
        @test sum(parse(Int,s["dimension"]) for s in sp["sectors"])==4^parse(Int,sp["L"])
        @test all(length(s["eigenvalues"])==parse(Int,s["dimension"]) for s in sp["sectors"])
    end
end

result_path=joinpath(root,get(ARGS,2,"results/audit_2026-10-11_check.json"))
@assert !isfile(result_path) "Choose a new checker-result filename"
write_result(result_path,(;status="passed",seconds=time()-started,
    test_counts=string(Test.get_test_counts(checks)),
    original_certificate_sha256=sha256_file(original),replay_sha256=sha256_file(p),
    source_hashes=source_hashes(root;dirs=("src",),files=["check_replay.jl","replay_audit.jl"]),
    julia=string(VERSION)))
