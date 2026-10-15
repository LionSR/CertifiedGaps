# NON-CERTIFYING Goerisch/Lehmann continuation experiment.
include("../src/common.jl")
using .AnchorStudy, LinearAlgebra, SparseArrays, Arpack, Random
using CertifiedGaps.Results
const A=AnchorStudy

function lowest_fixed(H,k)
    d=size(H,1)
    if d<=1300
        F=eigen(Symmetric(Matrix(H)),1:k)
        return F.values,F.vectors
    end
    # Asking for exactly the size of a degenerate multiplet can miss a direction
    # in scalar Lanczos. Oversample small spaces, as in the independent ring scan.
    nev=k<=8 ? min(d-2,20) : k
    Random.seed!(d+k)
    vals,V,nconv,_,_,_=eigs(H;nev,ncv=min(d,max(2nev+8,80)),which=:SR,
        tol=1e-10,maxiter=5000,v0=randn(d))
    @assert nconv==nev
    p=sortperm(vals)
    return vals[p[1:k]],V[:,p[1:k]]
end

function continuation(n,g; sector=0,mode="uniform")
    allp=read_exact_json(joinpath(A.RESULTS,"count_optimized.json"))["rows"]
    p=only(filter(x->A.ival(x["n"])==n,allp))
    u,v,c=A.num.((p["u"],p["v"],p["c"]))
    cut=A.tail_cutoff(n)
    if mode=="binned"
        bs=read_exact_json(joinpath(A.RESULTS,"binned_tails.json"))["rows"]
        cut=A.num(only(filter(x->A.ival(x["n"])==n,bs))["rho"])
    end
    sp=A.comparison_spectra([n÷2,n÷2],u,v;c)[sector]
    K=searchsortedfirst(sp,cut)-1
    initialK=K
    rho=sp[K+1]-2e-6
    @assert sp[K]<rho
    H0=A.comparison_matrix(n,sector,u,v,c)
    H=Float64.(A.ring(n,sector;twist=g=="I" ? 1 : -1))
    D=H-H0
    dnorm=2maximum(eigmax(Symmetric(X.H-u*X.L-v*X.Q))-c for X in A.LOCAL)
    t=0.; dt=min(0.1,0.8*(rho-sp[K])/dnorm)
    records=Any[]; attempts=0; started=time()
    suffix=mode=="uniform" ? "" : "_binned"
    path=joinpath(A.RESULTS,"homotopy_n$(n)_g$(g)_m$(sector)$suffix.json")
    if isfile(path) && get(ENV,"ANCHOR_FRESH","0")!="1"
        old=read_exact_json(path)
        if old["complete"]
            println("Already complete: $path")
            return
        end
        records=old["stages"]
        t=A.num(last(records)["t"]);K=A.ival(last(records)["k_out"])
        rho=A.num(last(records)["rho_out"]);attempts=A.ival(old["attempts"])
        started=time()-A.num(old["seconds"]);dt=0.1
        println("Resuming n=$n g=$g m=$sector at t=$t K=$K")
    end
    while t<1-1e-12 && K>0
        next=min(1.,t+dt)
        Ht=H0+next*D
        tick=time(); vals,V=lowest_fixed(Ht,K); attempts+=1
        margin=rho-last(vals)
        if margin<1e-8
            dt/=2
            @assert dt>1e-10 "homotopy stalled"
            continue
        end
        W=Ht*V-rho*V
        M=Symmetric(transpose(V)*W)
        N=Symmetric(transpose(W)*W)
        τ=eigvals(M,N)
        @assert maximum(τ)<0
        lower=sort(rho .+ 1 ./ τ) .- 1e-7
        # Keep the moving anchor well above the final target. Dropping at cut+epsilon
        # collapses the available margin and can stall at a degenerate crossing.
        drop_threshold=cut+0.9*(rho-cut)
        newK=searchsortedlast(lower,drop_threshold)
        next_rho=newK<K ? min(rho,lower[newK+1]-1e-7) : rho
        residual=maximum(norm(Ht*V[:,i]-vals[i]*V[:,i]) for i in 1:K)
        push!(records,Dict("stage"=>length(records)+1,"t"=>next,"k_in"=>K,
            "k_out"=>newK,"rho_in"=>rho,"rho_out"=>next_rho,
            "margin"=>margin,"residual"=>residual,"max_lower_error"=>
            maximum(vals-lower),"seconds"=>time()-tick))
        println("n=$n g=$g m=$sector stage=$(length(records)) t=$(round(next;digits=6)) "*
            "K=$K->$newK margin=$margin elapsed=$(round(time()-started;digits=1))")
        flush(stdout)
        K=newK;rho=next_rho;t=K==0 ? 1. : next
        # A numerical predictor only: every accepted stage still checks M<0.
        slope=K>0 ? dot(V[:,K],D*V[:,K]) : 1.
        dt=K>0 ? min(0.2,0.98*(rho-vals[K])/max(slope,1e-8)) : 1.
        write_result(path,Dict("metadata"=>A.metadata(Dict("n"=>n,"twist"=>g,
            "sector"=>sector,"mode"=>mode,
            "type"=>"Float64 generalized-eigenvalue continuation")),
            "n"=>n,"twist"=>g,"Sz"=>sector,"target"=>cut,"initial_k"=>initialK,
            "final_k"=>K,"stages"=>records,"attempts"=>attempts,
            "complete"=>t>=1-1e-12,"seconds"=>time()-started,"D_norm"=>dnorm))
        @assert length(records)<500
    end
end
if abspath(PROGRAM_FILE)==(@__FILE__)
    n=isempty(ARGS) ? 8 : parse(Int,ARGS[1])
    g=length(ARGS)>=2 ? ARGS[2] : "I"
    mode=length(ARGS)>=3 ? ARGS[3] : "uniform"
    sector=length(ARGS)>=4 ? parse(Int,ARGS[4]) : 0
    continuation(n,g;mode,sector)
end
