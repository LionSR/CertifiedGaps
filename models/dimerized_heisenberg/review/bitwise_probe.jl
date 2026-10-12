# Independent acceptance run. Never overwrites the production certificates.
include("../check_verified.jl")
const QQReview=Rational{BigInt}
const CQReview=Complex{QQReview}
@testset "Independent acceptance: full reruns and exact eight-spin residuals" begin
    check_verified_results()
    paths=sort(filter(p->startswith(basename(p),"verified_delta_")&&endswith(p,".json"),
        readdir(joinpath(@__DIR__,"../results");join=true)))
    for path in paths
        data=JSON3.read(read(path,String)); delta=parseq(data.delta); beta=parseq(data.beta0)
        println("ACCEPTANCE RERUN delta=",delta); flush(stdout)
        fresh=Dict{Int,Any}()
        for stored in data.rings
            L=Int(stored.spins)
            ring,E,r=verified_ring(delta,L)
            @test ring["dimension"]==stored.dimension
            @test length(ring["blocks"])==length(stored.blocks)
            for (a,b) in zip(ring["blocks"],stored.blocks)
                for key in ("dimension","down_spins","cell_momentum","spin_flip_parity",
                    "multiplicity","h_component_max","q_component_max","lambda_max_abs",
                    "raw_residual_component_max","raw_orthogonality_component_max",
                    "matrix_error_upper","residual_norm_upper","orthogonality_norm_upper",
                    "eigenvalue_radius","eigenvalues")
                    @test a[key]==b[Symbol(key)]
                end
            end
            fresh[L]=(E,r)
        end
        zs=enclosed_logZ(fresh[8]...,beta)
        zl=enclosed_logZ(fresh[16]...,beta)
        zd=enclosed_logZ(fresh[16]...,2beta)
        es=rational_upper(1-exp(zl-2zs)); et=rational_upper(1-exp(zd-2zl))
        @test es==parseq(data.epsilon_s_upper)
        @test et==parseq(data.epsilon_t_upper)
        @test certified_gap_bound(es,et,beta).bound==data.recurrence_bound
        # Independently recompute the residual and Gram defect with exact rational
        # matrix multiplication for every eight-spin sector, not the stored radii.
        exactblocks=0
        for q in 0:4
            orbits,flip=VerifiedSpectra.group_orbits(8,q)
            for k in 0:3, parity in (flip ? (1,-1) : (1,))
                built=VerifiedSpectra.verified_block(delta,8,q,k,parity,orbits)
                built===nothing && continue
                H,herr=built; d=size(H,1)
                rec,vals,radius=VerifiedSpectra.enclose_block(H,herr,8)
                eig=eigen(Hermitian(H)); Q=CQReview.(eig.vectors)
                residual=CQReview.(H)*Q-Q*Diagonal(QQReview.(eig.values))
                defect=Q'*Q-Matrix{CQReview}(I,d,d)
                @test sum(abs2,residual)<=parseq(rec["residual_norm_upper"])^2
                @test sum(abs2,defect)<=parseq(rec["orthogonality_norm_upper"])^2
                exactblocks+=1
            end
        end
        println("ACCEPTED delta=",delta," bound=",data.recurrence_bound,
            " exact-eight-spin-blocks=",exactblocks); flush(stdout)
    end
end
