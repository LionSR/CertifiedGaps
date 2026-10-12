using ThermalMPOGap, ThermalMPOGap.Results, IntervalArithmetic, LinearAlgebra
include("verified_spectra.jl")
include("rp.jl")
using .VerifiedSpectra, .DimerizedRP
BLAS.set_num_threads(1)
function run_verified(delta::Rational{BigInt},beta::Rational{BigInt})
    @assert 0<=delta<1 && beta>0
    a=audit_rp(delta); @assert a.reflected_real_identity
    t=time()
    small,Esmall,rsmall=verified_ring(delta,8)
    large,Elarge,rlarge=verified_ring(delta,16)
    lsmall=enclosed_logZ(Esmall,rsmall,beta)
    llarge=enclosed_logZ(Elarge,rlarge,beta)
    ldouble=enclosed_logZ(Elarge,rlarge,2beta)
    ps=exp(llarge-2lsmall); pt=exp(ldouble-2llarge)
    epss=rational_upper(1-ps); epst=rational_upper(1-pt)
    @assert 0<=epss<1 && 0<=epst<1
    bound=certified_gap_bound(epss,epst,beta)
    @assert bound.bound!==nothing
    root=normpath(joinpath(@__DIR__,"../.."))
    sources=source_hashes(root;dirs=("julia/src",),files=(
        "certificates/dimerized_heisenberg/rp.jl",
        "certificates/dimerized_heisenberg/numerical_ed.jl",
        "certificates/dimerized_heisenberg/verified_spectra.jl",
        "certificates/dimerized_heisenberg/run_verified.jl",
        "certificates/dimerized_heisenberg/check_verified.jl",
        "certificates/dimerized_heisenberg/test/verified_independent.jl",
        "certificates/dimerized_heisenberg/ENCLOSURE_PROOF.md",
        "julia/Project.toml","julia/Manifest.toml",
        "outputs/thermal-mpo-gap-note.tex","outputs/thermal-mpo-gap-scope.tex"))
    record=Dict("metadata"=>run_metadata(Dict("delta"=>string(delta),"beta0"=>string(beta),"n0_cells"=>4);
            source_sha256=sources,packages=(ThermalMPOGap,IntervalArithmetic),
            status="verified_spectral_seeds_general_reflected_RP",certified=true),
        "certification_scope"=>"General reflected RP theorem Theta=F_cell K; NOT the original stronger real identical-factor representation",
        "original_stronger_RP_condition_satisfied"=>false,
        "general_reflected_RP_exact"=>true,"delta"=>string(delta),"n0_cells"=>4,"beta0"=>string(beta),
        "epsilon_s_upper"=>string(epss),"epsilon_t_upper"=>string(epst),
        "log_partition_enclosures"=>Dict("n_beta"=>string.([rational_lower(lsmall),rational_upper(lsmall)]),
            "2n_beta"=>string.([rational_lower(llarge),rational_upper(llarge)]),
            "2n_2beta"=>string.([rational_lower(ldouble),rational_upper(ldouble)])),
        "recurrence_bound"=>bound.bound,"recurrence_bound_exact_binary"=>string(Rational{BigInt}(bound.bound)),
        "recurrence_diagnostics"=>bound,"validity_range"=>"even L >= 4 cells; 8,12,16,... physical spins",
        "seconds_total"=>time()-t,"rings"=>[small,large])
    path=joinpath(@__DIR__,"results","verified_delta_$(numerator(delta))_$(denominator(delta)).json")
    write_result(path,record)
    println("SAVED $path epsilon_s=$(Float64(epss)) epsilon_t=$(Float64(epst)) bound=$(bound.bound) seconds=$(time()-t)")
end
if abspath(PROGRAM_FILE)==@__FILE__
    delta=length(ARGS)>0 ? parse(BigInt,ARGS[1])//parse(BigInt,ARGS[2]) : big(1)//2
    beta=delta==1//2 ? big(6)//1 : delta==3//10 ? big(13)//2 : big(29)//4
    run_verified(delta,beta)
end
