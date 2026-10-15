# Checker for the dimerized spin-1 certificates.
#
#     julia --project=models/spin1_dimerized models/spin1_dimerized/check.jl          # ~1 min
#     julia --project=models/spin1_dimerized models/spin1_dimerized/check.jl --full   # + 17-36 min
#
# The fast check (test/runtests.jl) recomputes the seeds and bounds from the stored log Z
# enclosures and recertifies delta = 9/10 from scratch. `--full` also recertifies
# delta = 7/10 (twelve-spin ring) into a temporary file and checks that its log Z enclosures
# overlap the stored ones and that it proves a positive bound.
include("test/runtests.jl")
if "--full" in ARGS
    @testset "Fresh recertification of delta = 7/10" begin
        stored = read_exact_json(joinpath(ROOT, "results", "certificate_d7.json"))
        path = tempname() * ".json"
        fresh = certify_selected(7//10; n0 = 3, beta0 = 17//4, output_path = path)
        rm(path; force = true)
        @test fresh["bound"] > 0
        for key in ("logZ_n_beta", "logZ_2n_beta", "logZ_2n_2beta")
            a = parserat.((stored[key]["lower"], stored[key]["upper"]))
            b = parserat.((fresh[key]["lower"], fresh[key]["upper"]))
            @test max(a[1], b[1]) <= min(a[2], b[2])
        end
        @test abs(fresh["bound"] - stored["bound"]) < 1e-8
    end
end
