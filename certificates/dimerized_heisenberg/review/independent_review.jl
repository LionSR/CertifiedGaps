# Independent acceptance: fresh solves, rigorous refined radii, stored-interval containment.
# This does not require bitwise-identical eigenvectors and never overwrites production data.
include("../check_verified.jl")
const QQReview=Rational{BigInt}
const CQReview=Complex{QQReview}
function refined_witness(rec,L)
    d=Int(rec["dimension"])
    h=parseq(rec["matrix_error_upper"])
    q=parseq(rec["q_component_max"])
    lam=parseq(rec["lambda_max_abs"])
    rawr=parseq(rec["raw_residual_component_max"])
    rawe=parseq(rec["raw_orthogonality_component_max"])
    coarse=VerifiedSpectra.scalar_witness(d,L,parseq(rec["h_component_max"]),
        q,lam,rawr,rawe,h)
    @assert coarse.ebound<1 # Thus ||Q||<2, independently of the refinement.
    a=big(3)*L//4
    # Cauchy on two concatenated real dot-product vectors gives these sums.
    rerr=min(coarse.rerr,coarse.gamma*(2*(a+h)+q*lam)+coarse.underflow)
    eerr=min(coarse.eerr,coarse.gamma*5+coarse.underflow)
    # ||complex dxd matrix||F <= sqrt(2)*d*componentmax < (3/2)*d*componentmax.
    # This also strictly improves the coarse radius of scalar/small blocks.
    factor=big(3)*d//2
    r=factor*(rawr+rerr); e=factor*(rawe+eerr)
    @assert e<1
    return (;r,e,radius=r+2*h+(a+lam)*e)
end
reviewrecords=Any[]
elapsed=@elapsed @testset "Independent acceptance: containment from fresh full spectra" begin
    check_verified_results()
    paths=sort(filter(p->startswith(basename(p),"verified_delta_")&&endswith(p,".json"),
        readdir(joinpath(@__DIR__,"../results");join=true)))
    for path in paths
        data=JSON3.read(read(path,String)); delta=parseq(data.delta); beta=parseq(data.beta0)
        println("ACCEPTANCE RERUN delta=",delta); flush(stdout)
        fresh=Dict{Int,Any}(); rings=Any[]; mismatches=Dict{String,Int}()
        for stored in data.rings
            L=Int(stored.spins)
            ring,E,r=verified_ring(delta,L)
            @test ring["dimension"]==stored.dimension
            @test length(ring["blocks"])==length(stored.blocks)
            centerlist=Float64[]; refinedlist=QQReview[]; blocks=Any[]
            for (a,b) in zip(ring["blocks"],stored.blocks)
                for key in ("dimension","down_spins","cell_momentum","spin_flip_parity","multiplicity")
                    @test a[key]==b[Symbol(key)]
                end
                for key in ("h_component_max","q_component_max","lambda_max_abs",
                    "raw_residual_component_max","raw_orthogonality_component_max",
                    "matrix_error_upper","residual_norm_upper","orthogonality_norm_upper",
                    "eigenvalue_radius","eigenvalues")
                    a[key]==b[Symbol(key)] || (mismatches[key]=get(mismatches,key,0)+1)
                end
                w=refined_witness(a,L)
                av=parseq.(a["eigenvalues"]); bv=parseq.(b.eigenvalues)
                oldradius=parseq(b.eigenvalue_radius)
                displacement=maximum(abs.(av.-bv))
                margin=oldradius-displacement-w.radius
                # Ordered Weyl intervals independently regenerated around fresh centers
                # lie wholly inside the production intervals, including multiplicities.
                @test margin>=0
                for _ in 1:Int(a["multiplicity"])
                    append!(centerlist,Float64.(av)); append!(refinedlist,fill(w.radius,length(av)))
                end
                push!(blocks,Dict("down_spins"=>a["down_spins"],"cell_momentum"=>a["cell_momentum"],
                    "spin_flip_parity"=>a["spin_flip_parity"],"dimension"=>a["dimension"],
                    "multiplicity"=>a["multiplicity"],"fresh_witness"=>a,
                    "refined_residual_norm"=>string(w.r),"refined_orthogonality_norm"=>string(w.e),
                    "refined_eigenvalue_radius"=>string(w.radius),
                    "maximum_center_displacement"=>string(displacement),
                    "minimum_stored_interval_containment_margin"=>string(margin)))
            end
            @test length(centerlist)==2^L
            fresh[L]=(centerlist,refinedlist)
            push!(rings,Dict("spins"=>L,"dimension"=>ring["dimension"],"blocks"=>blocks,
                "blas_config"=>ring["blas_config"],"seconds"=>ring["seconds"]))
        end
        zs=enclosed_logZ(fresh[8]...,beta); zl=enclosed_logZ(fresh[16]...,beta)
        zd=enclosed_logZ(fresh[16]...,2beta)
        es=rational_upper(1-exp(zl-2zs)); et=rational_upper(1-exp(zd-2zl))
        @test es<=parseq(data.epsilon_s_upper)
        @test et<=parseq(data.epsilon_t_upper)
        freshbound=certified_gap_bound(es,et,beta).bound
        @test freshbound!==nothing && freshbound>=data.recurrence_bound
        exactblocks=0
        for q in 0:4
            orbits,flip=VerifiedSpectra.group_orbits(8,q)
            for k in 0:3, parity in (flip ? (1,-1) : (1,))
                built=VerifiedSpectra.verified_block(delta,8,q,k,parity,orbits)
                built===nothing && continue
                H,herr=built; d=size(H,1)
                rec,vals,radius=VerifiedSpectra.enclose_block(H,herr,8)
                w=refined_witness(rec,8)
                eig=eigen(Hermitian(H)); Q=CQReview.(eig.vectors)
                residual=CQReview.(H)*Q-Q*Diagonal(QQReview.(eig.values))
                defect=Q'*Q-Matrix{CQReview}(I,d,d)
                @test sum(abs2,residual)<=w.r^2
                @test sum(abs2,defect)<=w.e^2
                exactblocks+=1
            end
        end
        push!(reviewrecords,Dict("delta"=>string(delta),"production_certificate"=>basename(path),
            "production_sha256"=>sha256_file(path),"beta"=>string(beta),
            "stored_epsilon_s"=>String(data.epsilon_s_upper),"stored_epsilon_t"=>String(data.epsilon_t_upper),
            "fresh_refined_epsilon_s"=>string(es),"fresh_refined_epsilon_t"=>string(et),
            "accepted_stored_bound"=>String(data.recurrence_bound_exact_binary),
            "fresh_refined_bound"=>string(QQReview(freshbound)),
            "non_bitwise_equal_block_counts"=>mismatches,"exact_rational_eight_spin_blocks"=>exactblocks,
            "all_stored_eigenvalue_intervals_independently_contained"=>true,"rings"=>rings))
        println("ACCEPTED delta=",delta," stored_bound=",data.recurrence_bound,
            " fresh_refined_bound=",freshbound," exact-eight-spin-blocks=",exactblocks,
            " mismatch_counts=",mismatches); flush(stdout)
    end
end
root=normpath(joinpath(@__DIR__,"../../.."))
files=vcat(["certificates/dimerized_heisenberg/review/independent_review.jl",
    "certificates/dimerized_heisenberg/verified_spectra.jl",
    "certificates/dimerized_heisenberg/ENCLOSURE_PROOF.md"],
    ["certificates/dimerized_heisenberg/results/"*r["production_certificate"] for r in reviewrecords])
write_result(joinpath(@__DIR__,"acceptance.json"),Dict(
    "metadata"=>run_metadata(Dict("scope"=>"independent_general_reflected_RP_acceptance");
        source_sha256=source_hashes(root;dirs=(),files=Tuple(files)),
        packages=(ThermalMPOGap,IntervalArithmetic),status="pass",certified=true),
    "original_strict_real_identical_hypothesis"=>"FAIL",
    "general_reflected_RP_certificates"=>"PASS",
    "backend_assumption"=>"Classical binary64 CPU OpenBLAS GEMM with ordinary multiply/add/FMA",
    "seconds"=>elapsed,"certificates"=>reviewrecords))
