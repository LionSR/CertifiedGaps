# All spectra and derived thermal values are NON-CERTIFYING Float64 calculations.
include(joinpath(@__DIR__, "..", "src", "common.jl"))
using SparseArrays, KrylovKit, Random, Serialization
BLAS.set_num_threads(1)

"Spin-1 Sz block with a physical rotation at the closing bond."
function hamiltonian(n,m,twist)
    strides=3 .^ (0:n-1)
    states=[s for s in 0:3^n-1 if sum(digits(s;base=3,pad=n))==n+m]
    lookup=Dict(s=>i for (i,s) in enumerate(states))
    rows=Int[]; cols=Int[]
    T=twist==2 ? ComplexF64 : Float64
    vals=T[]
    phase=twist==2 ? cis(2pi/3) : twist==1 ? -1.0 : 1.0
    for (i,s) in enumerate(states)
        spins=digits(s;base=3,pad=n).-1
        diag=0
        for j in 1:n
            k=mod1(j+1,n)
            diag+=spins[j]*spins[k]
            for delta in (-1,1)
                abs(spins[j]-delta)<=1 && abs(spins[k]+delta)<=1 || continue
                target=s+delta*(strides[k]-strides[j])
                push!(rows,lookup[target]); push!(cols,i)
                push!(vals,k==1 ? phase^delta : one(T))
            end
        end
        push!(rows,i);push!(cols,i);push!(vals,diag)
    end
    H=sparse(rows,cols,vals,length(states),length(states))
    @assert norm(H-H')<1e-10
    H
end

function spectrum(n,twist; nev=96)
    sectors=Any[]
    for m in 0:n
        H=hamiltonian(n,m,twist)
        dim=size(H,1); started=time()
        if dim<=1400
            es=eigvals(Hermitian(Matrix(H)))
            method="dense Float64 ED"
            residual=0.0
        else
            Random.seed!(10000n+100twist+m)
            take=min(nev,dim-2)
            es,vs,info=eigsolve(H,randn(eltype(H),dim),take,:SR;
                ishermitian=true,krylovdim=min(dim,4take+40),tol=1e-11,maxiter=200)
            info.converged>=take || error("unconverged Lanczos n=$n g=$twist m=$m")
            es=real.(es[1:take]); residual=maximum(info.normres[1:take])
            method="Float64 Lanczos; no certified counting anchor"
        end
        push!(sectors,Dict("m"=>m,"dim"=>dim,"eigenvalues"=>es,
            "multiplicity"=>m==0 ? 1 : 2,"method"=>method,
            "maximum_ritz_residual"=>residual,"seconds"=>time()-started))
        println("spectrum n=$n g=$twist m=$m dim=$dim count=$(length(es)) ",
            "seconds=$(round(time()-started;digits=2))")
        flush(stdout)
    end
    @assert sum(s["dim"]*s["multiplicity"] for s in sectors)==3^n
    Dict("n"=>n,"twist"=>twist,"certified"=>false,"sectors"=>sectors)
end

function main()
    mkpath(CACHE)
    lengths=isempty(ARGS) ? (4:10) : parse.(Int,ARGS)
    for n in lengths, g in 0:2
        n == 12 && g == 2 && continue  # Baseline workload needs only the two real twists.
        path=joinpath(CACHE,"spectrum-$n-$g.json")
        isfile(path) && continue
        result=spectrum(n,g;nev=n==12 ? 48 : 96)
        result["metadata"]=metadata("Float64 spectra; n=$n g=$g")
        write_compact(path,result)
    end
end
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
