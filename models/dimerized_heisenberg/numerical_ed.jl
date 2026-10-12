# Full dense symmetry-block spectra. These Float64 eigenvalues are NOT certified.
using ThermalMPOGap, LinearAlgebra
using ThermalMPOGap.Results, ThermalMPOGap.ExactDiagonalization, ThermalMPOGap.Seeds
BLAS.set_num_threads(1)
const ROOT_ED = normpath(joinpath(@__DIR__, "../.."))
rotate_cell(s,L) = ((s << 2) & ((1<<L)-1)) | (s >> (L-2))

function group_orbits(L, q)
    n=L÷2; mask=(1<<L)-1
    states=[s for s in 0:mask if count_ones(s)==q]
    seen=falses(1<<L); orbits=Any[]
    useflip=(2q==L)
    for rep in states
        seen[rep+1] && continue
        images=Dict{Int,Tuple{Int,Int}}(); stabilizers=Tuple{Int,Int}[]
        for f in 0:Int(useflip)
            s=f==0 ? rep : xor(rep,mask)
            for t in 0:n-1
                if haskey(images,s)
                    t0,f0=images[s]
                    push!(stabilizers,(t-t0,f-f0))
                else
                    images[s]=(t,f); seen[s+1]=true
                end
                s=rotate_cell(s,L)
            end
        end
        push!(orbits,(;rep,images,stabilizers))
    end
    @assert sum(length(o.images) for o in orbits)==binomial(L,q)
    return orbits,useflip
end

function dense_spectrum(delta,L)
    n=L÷2; levels=Float64[]; records=Any[]; dimension_sum=0
    for q in 0:L÷2
        orbits,useflip=group_orbits(L,q)
        mult=useflip ? 1 : 2
        qdim=0
        for k in 0:n-1, parity in (useflip ? (1,-1) : (1,))
            odd=parity==-1 ? 1 : 0
            allowed=[o for o in orbits if all(mod(2*k*t+n*odd*f,2*n)==0 for (t,f) in o.stabilizers)]
            dim=length(allowed); qdim+=dim
            dim==0 && continue
            labels=zeros(Int,1<<L); coeff=zeros(ComplexF64,1<<L)
            for (i,o) in enumerate(allowed), (s,(t,f)) in o.images
                labels[s+1]=i
                coeff[s+1]=cispi(2*k*t/n)*parity^f/sqrt(length(o.images))
            end
            H=zeros(ComplexF64,dim,dim)
            for (j,o) in enumerate(allowed)
                s=o.rep; scale=sqrt(length(o.images))
                diag=0.0
                for site in 0:L-1
                    next=mod(site+1,L)
                    J=1+(iseven(site) ? delta : -delta)
                    opposite=((s>>site)&1)!=((s>>next)&1)
                    diag+=J*(opposite ? -0.25 : 0.25)
                    if opposite
                        target=xor(s,(1<<site)|(1<<next)); i=labels[target+1]
                        i==0 || (H[i,j]+=J/2*scale*conj(coeff[target+1]))
                    end
                end
                H[j,j]+=diag
            end
            defect=maximum(abs,H-H')
            @assert defect < 1e-10
            duration=@elapsed E=eigvals(Hermitian(H)) # Dense LAPACK only.
            for _ in 1:mult; append!(levels,E); end
            dimension_sum+=mult*dim
            push!(records,Dict("down_spins"=>q,"cell_momentum"=>k,
                "spin_flip_parity"=>(useflip ? parity : nothing),"dimension"=>dim,
                "multiplicity"=>mult,"hermiticity_defect"=>defect,"seconds"=>duration))
        end
        @assert qdim==binomial(L,q)
    end
    @assert dimension_sum==2^L
    @assert length(levels)==2^L
    return PartialSpectrum(sort!(levels),2^L,true),records
end

function seedrow(r)
    r===nothing && return nothing
    Dict("n0_cells"=>r.n₀,"beta0"=>r.β₀,"epsilon_s"=>r.εₛ,
        "epsilon_t"=>r.εₜ,"numerical_bound"=>r.bound)
end
function run_numerics(; maxseed=4)
    # n0=3 is retained as a diagnostic but excluded from theorem-eligible best seeds.
    ns=collect(2:maxseed); lengths=sort!(unique(vcat(ns,2 .* ns)))
    deltas=[big(1)//2,big(3)//10,big(1)//4,big(1)//5,big(3)//20,big(1)//10]
    for d in deltas
        delta=Float64(d); specs=Dict{Int,PartialSpectrum}(); sizes=Any[]
        total=@elapsed begin
            for cells in lengths
                elapsed=@elapsed spec,blocks=dense_spectrum(delta,2*cells)
                specs[cells]=spec
                push!(sizes,Dict("cells"=>cells,"spins"=>2*cells,"dimension"=>spec.N,
                    "complete_dense_spectrum"=>true,"ground_energy"=>spec.E[1],
                    "reference_gap"=>spec.E[2]-spec.E[1],"seconds"=>elapsed,"sectors"=>blocks))
                println("delta=$d cells=$cells dimension=$(spec.N) gap=$(spec.E[2]-spec.E[1]) seconds=$elapsed"); flush(stdout)
            end
        end
        grid=scan_seeds(spectral_defects(specs),ns,0.25:0.25:32.0)
        diagnostics=Any[]
        for n in ns
            rows=filter(r->r.n₀==n,grid)
            peak=argmax(r->min(1-r.εₛ,1-r.εₜ),rows)
            push!(diagnostics,Dict("n0_cells"=>n,"even_seed_eligible"=>iseven(n),
                "max_min_purity"=>min(1-peak.εₛ,1-peak.εₜ),"peak_beta"=>peak.β₀,
                "energy_difference"=>specs[2n].E[1]-2*specs[n].E[1],
                "best_seed"=>seedrow(best_seed(rows))))
        end
        best=best_seed(filter(r->iseven(r.n₀),grid))
        smallest=minimum((r.n₀ for r in grid if iseven(r.n₀)&&r.bound!==nothing);init=typemax(Int))
        refgap=specs[maximum(lengths)].E[2]-specs[maximum(lengths)].E[1]
        data=Dict("metadata"=>run_metadata(Dict("delta"=>string(d),"maxseed"=>maxseed,
                "beta_grid"=>"0.25:0.25:32.0");
                source_sha256=source_hashes(ROOT_ED;dirs=("julia/src",),files=("certificates/dimerized_heisenberg/numerical_ed.jl",)),
                packages=(ThermalMPOGap,),status="numerical_not_certified",certified=false),
            "delta"=>string(d),"seconds_spectra"=>total,"sizes"=>sizes,
            "diagnostics"=>diagnostics,"best_even_seed"=>seedrow(best),
            "smallest_successful_even_n0"=>(smallest==typemax(Int) ? nothing : smallest),
            "reference_gap"=>refgap,"reference_spins"=>2*maximum(lengths),
            "bound_to_finite_ring_gap_ratio"=>(best===nothing ? nothing : best.bound/refgap),
            "scan"=>[seedrow(r) for r in grid],
            "rigor_limit"=>"Full dense spectra with dimension checks, but no eigenvalue enclosures or certified purity seeds; all gap values numerical.")
        path=joinpath(@__DIR__,"results","numerical_delta_$(numerator(d))_$(denominator(d)).json")
        write_result(path,data)
        println("saved $path best=$(seedrow(best))"); flush(stdout)
    end
end
if abspath(PROGRAM_FILE)==@__FILE__
    run_numerics()
end
