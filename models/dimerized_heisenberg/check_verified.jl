using CertifiedGaps, CertifiedGaps.Results, JSON3, Test, IntervalArithmetic, LinearAlgebra
isdefined(@__MODULE__,:check_historical_hashes) || include("provenance.jl")
isdefined(@__MODULE__,:VerifiedSpectra) || include("verified_spectra.jl")
isdefined(@__MODULE__,:DimerizedRP) || include("rp.jl")
using .VerifiedSpectra, .DimerizedRP
BLAS.set_num_threads(1)
parseq(s)=(p=split(String(s),"//"); parse(BigInt,p[1])//parse(BigInt,p[2]))
function check_verified_file(path)
    data=JSON3.read(read_json_text(path))
    @test data.original_stronger_RP_condition_satisfied===false
    @test data.general_reflected_RP_exact===true
    @test 0<=parseq(data.delta)<1
    @test audit_rp(parseq(data.delta)).reflected_real_identity
    check_source_hashes(data.metadata)
    spectra=Dict{Int,Tuple{Vector{Float64},Vector{Rational{BigInt}}}}()
    for ring in data.rings
        E=Float64[]; radii=Rational{BigInt}[]; totaldim=0
        for block in ring.blocks
            d=Int(block.dimension)
            w=VerifiedSpectra.scalar_witness(d,Int(ring.spins),parseq(block.h_component_max),
                parseq(block.q_component_max),parseq(block.lambda_max_abs),
                parseq(block.raw_residual_component_max),parseq(block.raw_orthogonality_component_max),
                parseq(block.matrix_error_upper))
            @test w.steps*VerifiedSpectra.UNIT<1
            @test w.ebound==parseq(block.orthogonality_norm_upper)<1
            @test w.rbound==parseq(block.residual_norm_upper)
            @test w.radius==parseq(block.eigenvalue_radius)>0
            vals=parseq.(block.eigenvalues)
            @test length(vals)==d
            @test all(x->Rational{BigInt}(Float64(x))==x,vals)
            @test issorted(vals)
            for _ in 1:Int(block.multiplicity)
                append!(E,Float64.(vals)); append!(radii,fill(w.radius,d))
            end
            totaldim+=d*Int(block.multiplicity)
        end
        @test totaldim==2^Int(ring.spins)==Int(ring.dimension)==length(E)
        for q in 0:Int(ring.spins)÷2
            @test sum(Int(b.dimension) for b in ring.blocks if Int(b.down_spins)==q)==binomial(Int(ring.spins),q)
        end
        spectra[Int(ring.spins)]=(E,radii)
    end
    beta=parseq(data.beta0)
    small=enclosed_logZ(spectra[8]...,beta)
    large=enclosed_logZ(spectra[16]...,beta)
    double=enclosed_logZ(spectra[16]...,2beta)
    for (key,value) in ((:n_beta,small),(:var"2n_beta",large),(:var"2n_2beta",double))
        stored=parseq.(data.log_partition_enclosures[key])
        @test stored==[rational_lower(value),rational_upper(value)]
    end
    es=rational_upper(1-exp(large-2small)); et=rational_upper(1-exp(double-2large))
    @test es==parseq(data.epsilon_s_upper)
    @test et==parseq(data.epsilon_t_upper)
    result=certified_gap_bound(es,et,beta)
    @test result.bound!==nothing
    @test Rational{BigInt}(result.bound)==parseq(data.recurrence_bound_exact_binary)
    @test result.bound==data.recurrence_bound
    println("verified stored delta=$(data.delta) bound=$(result.bound)")
end
function check_verified_results()
    paths=filter(p->startswith(basename(p),"verified_delta_")&&endswith(p,".json.gz"),
        readdir(joinpath(@__DIR__,"results");join=true))
    @test length(paths)==3
    for path in sort(paths); check_verified_file(path); end
    check_certificate_summary(paths)
end
"The plain results/certificate.json must copy the seeds of the compressed certificates exactly."
function check_certificate_summary(paths)
    summary=JSON3.read(read(joinpath(@__DIR__,"results","certificate.json"),String))
    rows=Dict(String(r.source_file)=>r for r in summary.certificates)
    @test length(rows)==length(paths)
    for path in paths
        r=rows["results/"*basename(path)]; data=JSON3.read(read_json_text(path))
        @test String(r.source_file_sha256)==sha256_file(path)
        @test Int(r.n0_cells)==Int(data.n0_cells)==4 && iseven(Int(r.n0_cells))
        for key in (:delta,:beta0,:epsilon_s_upper,:epsilon_t_upper,:recurrence_bound_exact_binary)
            @test String(r[key])==String(data[key])
        end
    end
end
if abspath(PROGRAM_FILE)==@__FILE__
    @testset "Verified spectral witnesses, partition sums and recurrence" begin
        check_verified_results()
    end
end
