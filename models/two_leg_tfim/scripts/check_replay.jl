# Checks a replay written by replay_audit.jl against the stored certificate, recomputing the
# comparison claims, log Z, the rational seeds and the bound rather than trusting saved flags.
#     julia --project=models/two_leg_tfim models/two_leg_tfim/scripts/check_replay.jl \
#         results/replay.json results/replay_check.json
using Test, SHA, CertifiedGaps, CertifiedGaps.Results
const ROOT = dirname(@__DIR__)
include(joinpath(ROOT, "src", "Ladder.jl"))
using .Ladder
include(joinpath(ROOT, "src", "Provenance.jl"))
using .Provenance
root=ROOT
length(ARGS)>=1 || error("usage: check_replay.jl <replay.json> [<check-output.json>]")
p=joinpath(root,ARGS[1])
c=read_exact_json(p)
started=time()
original=joinpath(root,ORIGINAL_CERTIFICATE)
old=read_certificate(original)
checks=@testset "Fresh spectral replay artifact" begin
    @test bytes2hex(sha256(certificate_bytes(original)))==c["original_certificate_sha256"]
    # A replay records the hashes of the sources it ran with; compare with current files.
    for (path,hash) in c["source_hashes"]
        @test source_sha256(root,path)==hash
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
    specs=stored_spectra(c)
    (;es,et)=seed_defects(specs[1],specs[2],beta)
    @test es==parse(Rational{BigInt},c["epsilon_s"])
    @test et==parse(Rational{BigInt},c["epsilon_t"])
    @test certified_gap_bound(es,et,beta).bound==c["bound"]
    @test Rational{BigInt}(c["bound"])==parse(Rational{BigInt},c["bound_exact_binary"])
    for sp in c["spectra"]
        @test sum(parse(Int,s["dimension"]) for s in sp["sectors"])==4^parse(Int,sp["L"])
        @test all(length(s["eigenvalues"])==parse(Int,s["dimension"]) for s in sp["sectors"])
    end
end

result_path=joinpath(root,get(ARGS,2,"results/replay_check.json"))
@assert !isfile(result_path) "Choose a new checker-result filename"
write_result(result_path,(;status="passed",seconds=time()-started,
    test_counts=string(Test.get_test_counts(checks)),
    original_certificate_sha256=bytes2hex(sha256(certificate_bytes(original))),
    replay_sha256=sha256_file(p),
    source_hashes=source_hashes(root;dirs=("src","scripts")),
    julia=string(VERSION)))
