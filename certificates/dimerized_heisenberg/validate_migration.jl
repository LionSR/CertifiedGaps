# Separate current-core validation; never rewrites historical results or source hashes.
using CertifiedGaps, CertifiedGaps.Results, JSON3, Test, SHA, LinearAlgebra
using IntervalArithmetic
BLAS.set_num_threads(1)
const HERE = @__DIR__
const ROOT = normpath(joinpath(HERE, "../.."))
const PIN = "ef0a3e6cf54ec56dcef5e8a37963d0184f6fd31b"
const CACHE = joinpath(HERE, ".compat")
mkpath(CACHE)
# The only migration applied to these two source copies is the package name.
for name in ("verified_spectra.jl", "numerical_ed.jl")
    write(joinpath(CACHE,name),replace(read(joinpath(HERE,name),String),
        "ThermalMPOGap"=>"CertifiedGaps"))
end
include(joinpath(CACHE,"verified_spectra.jl"))
using .VerifiedSpectra
include("rp.jl")
q(s) = (p=split(String(s),"//"); parse(BigInt,p[1])//parse(BigInt,p[2]))
records = Any[]
t0 = time()
@testset "Current CertifiedGaps core with immutable historical spectra" begin
    for name in ("verified_delta_1_2.json","verified_delta_3_10.json","verified_delta_1_4.json")
        data=JSON3.read(read(joinpath(HERE,"results",name),String))
        # Verify original provenance against the original git objects, not new core files.
        for (file,digest) in data.metadata.source_sha256
            path=String(file)
            bytes=startswith(path,"certificates/") ? read(joinpath(ROOT,path)) :
                read(`git -C $ROOT show $(PIN*":"*path)`)
            @test bytes2hex(sha256(bytes))==digest
        end
        @test DimerizedRP.audit_rp(q(data.delta)).reflected_real_identity
        rings=Dict()
        for ring in data.rings
            centers=Float64[]; radii=Rational{BigInt}[]
            for b in ring.blocks
                w=scalar_witness(Int(b.dimension),Int(ring.spins),q(b.h_component_max),
                    q(b.q_component_max),q(b.lambda_max_abs),q(b.raw_residual_component_max),
                    q(b.raw_orthogonality_component_max),q(b.matrix_error_upper))
                @test w.radius==q(b.eigenvalue_radius)
                @test w.ebound==q(b.orthogonality_norm_upper)<1
                @test w.rbound==q(b.residual_norm_upper)
                vals=q.(b.eigenvalues)
                @test all(x->Rational{BigInt}(Float64(x))==x, vals)
                @test length(vals)==Int(b.dimension)
                for _ in 1:Int(b.multiplicity)
                    append!(centers,Float64.(vals)); append!(radii,fill(w.radius,length(vals)))
                end
            end
            @test length(centers)==2^Int(ring.spins)==Int(ring.dimension)
            rings[Int(ring.spins)]=(centers,radii)
        end
        beta=q(data.beta0)
        a=enclosed_logZ(rings[8]...,beta)
        b=enclosed_logZ(rings[16]...,beta)
        c=enclosed_logZ(rings[16]...,2beta)
        es=rational_upper(1-exp(b-2a)); et=rational_upper(1-exp(c-2b))
        @test es==q(data.epsilon_s_upper)
        @test et==q(data.epsilon_t_upper)
        bound=certified_gap_bound(es,et,beta).bound
        @test Rational{BigInt}(bound)==q(data.recurrence_bound_exact_binary)
        push!(records,Dict("delta"=>String(data.delta),"epsilon_s"=>string(es),
            "epsilon_t"=>string(et),"bound_exact_binary"=>string(Rational{BigInt}(bound)),
            "historical_result_sha256"=>sha256_file(joinpath(HERE,"results",name))))
    end
end
hashes=source_hashes(ROOT;dirs=("julia/src",),files=("julia/Project.toml","julia/Manifest.toml",
    "certificates/dimerized_heisenberg/validate_migration.jl",
    "certificates/dimerized_heisenberg/verified_spectra.jl",
    "certificates/dimerized_heisenberg/numerical_ed.jl",
    "certificates/dimerized_heisenberg/rp.jl"))
write_result(joinpath(HERE,"results","migration_validation.json"),Dict(
    "metadata"=>run_metadata(Dict("historical_core_commit"=>PIN);
        source_sha256=hashes,packages=(CertifiedGaps,IntervalArithmetic),
        status="current_core_reconstruction_pass",certified=true),
    "validated_core_tree"=>strip(read(`git -C $ROOT rev-parse HEAD:julia`,String)),
    "validated_checkout_head"=>strip(read(`git -C $ROOT rev-parse HEAD`,String)),
    "transformation"=>"Replace ThermalMPOGap with CertifiedGaps in ignored copies of verified_spectra.jl and numerical_ed.jl only",
    "scope"=>"Exact RP, historical provenance, scalar witnesses, complete partition sums, rational seeds and recurrence; no new dense eigensolves",
    "seconds"=>time()-t0,"certificates"=>records))
