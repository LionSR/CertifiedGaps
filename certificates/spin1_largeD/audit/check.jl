# Fast checker of separately recorded fresh-replay evidence; does not redo dense ED.
include("assemble_replay.jl")
using Test
@testset "Fresh spectral replay evidence and rational certificates" begin
    for (D,beta,published) in ((4//1,9//4,694694192//1000000000),
                              (1603//500,381//125,2213456//1000000000))
        tag="D$(numerator(D))_$(denominator(D))"
        stored=readjson(joinpath(auditroot,"results/certificate_$tag.json"))
        for (path,hash) in stored["source_hashes"]
            @test bytes2hex(sha256(read(joinpath(auditroot,path))))==hash
        end
        for (path,hash) in stored["input_hashes"]
            @test bytes2hex(sha256(read(joinpath(auditroot,path))))==hash
        end
        fresh=assemble(D,beta,published;save=false)
        for key in ("D","n0","beta","eps_s","eps_t","bound","published_bound_rational",
                    "j_admissible","logZ","finite_gap_intervals","dimension_L6","dimension_L12",
                    "sector_count_L12","comparison","original_rational_seeds_equal",
                    "core_hashes","core_tree","worker_seconds","input_hashes","source_hashes")
            @test fresh[key]==stored[key]
        end
        @test Rational{BigInt}(fresh["bound"])>=published
    end
end
println("Fast audit check trusts freshly recorded radii; replay_worker.jl recomputes them from model inputs.")
