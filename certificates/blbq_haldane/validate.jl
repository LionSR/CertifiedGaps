using CertifiedGaps, CertifiedGaps.Results, LinearAlgebra
BLAS.set_num_threads(1)
start=time()
include("test/runtests.jl")
seconds=time()-start
@assert seconds<120 "Fast checker exceeded the requested two-minute limit"
write_result(joinpath(@__DIR__,"results/validation_current_main.json"),(;
    metadata=run_metadata((;test="test/runtests.jl");
        source_sha256=source_hashes(@__DIR__;dirs=("test",),
            files=("validate.jl","model.jl","sectors.jl","su2.jl","dimensions.jl","result_io.jl")),
        status="tests passed; no uniform gap certificate"),seconds,
    checks=["exact RP polynomial coefficient identities","exact span determinant",
        "small complete spectra agree with unreduced ED","stored rational recurrence rejects candidates",
        "all dense sector multiplicities sum to 3^L","stored input hashes agree",
        "full-spectrum first and second moments match analytic infinite-temperature traces",
        "exact small highest-weight projector identities and complete SU2 spectra through L7",
        "independent integer character dimensions for highest-weight blocks",
        "empty momentum and zero-rank highest-weight sectors",
        "complete L14 ground multiplet and low-energy reference agreement when available"],
    uniform_gap_certified=false))
println("Validation completed in $seconds seconds")
