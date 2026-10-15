# Exact rational residuals and the analytic four-spin spectrum versus the block enclosures.
isdefined(@__MODULE__, :check_verified_results) || include(joinpath(@__DIR__, "..", "check.jl"))
@testset "Independent exact residual and analytic four-spin validation" begin
    QQ=Rational{BigInt}; CQ=Complex{QQ}
    for delta in (big(1)//2,big(3)//10,big(1)//4)
        allvals=Float64[]; allr=QQ[]
        for q in 0:2
            orbits,flip=VerifiedSpectra.group_orbits(4,q)
            for k in 0:1, parity in (flip ? (1,-1) : (1,))
                built=VerifiedSpectra.verified_block(delta,4,q,k,parity,orbits)
                built===nothing && continue
                H,herr=built; d=size(H,1)
                record,vals,radius=VerifiedSpectra.enclose_block(H,herr,4)
                decomposition=eigen(Hermitian(H))
                Q=CQ.(decomposition.vectors); Hq=CQ.(H); V=QQ.(decomposition.values)
                R=Hq*Q-Q*Diagonal(V)
                defect=Q'*Q-Matrix{CQ}(I,d,d)
                @test sum(abs2,R)<=parseq(record["residual_norm_upper"])^2
                @test sum(abs2,defect)<=parseq(record["orthogonality_norm_upper"])^2
                # Deliberately wrong zero bound must be rejected for nonzero residual.
                @test iszero(sum(abs2,R)) || !(sum(abs2,R)<=0)
                for _ in 1:(flip ? 1 : 2)
                    append!(allvals,vals); append!(allr,fill(radius,d))
                end
            end
        end
        @test length(allvals)==16
        order=sortperm(allvals); vals=allvals[order]; rad=allr[order]
        # Complete exact four-spin spectrum: two singlets, three triplets, one quintet.
        # Use 256-bit intervals, independently of the 53-bit block enclosures.
        setprecision(BigFloat,256) do
            x=interval(BigFloat,numerator(delta))/interval(BigFloat,denominator(delta))
            root=sqrt(1+3x^2)
            analytic=vcat([-1-root,-1+root],fill(interval(BigFloat,-1),3),
                fill(-x,3),fill(x,3),fill(interval(BigFloat,1),5))
            sort!(analytic;by=inf)
            for (center,radius,exact) in zip(vals,rad,analytic)
                err=max(abs(inf(exact)-BigFloat(center)),abs(sup(exact)-BigFloat(center)))
                @test err<BigFloat(radius)
            end
            beta=big(6)//1
            computed=enclosed_logZ(allvals,allr,beta)
            direct=log(sum(exp.(-6 .* analytic)))
            @test BigFloat(inf(computed))<=inf(direct)<=sup(direct)<=BigFloat(sup(computed))
        end
    end
end
