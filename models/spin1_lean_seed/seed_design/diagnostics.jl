include("search_new.jl")

function direct_seed(input,N)
    pwr=N/input.nmax
    mult=input.multiplicities
    tail(k)=sum(mult[j]*(input.masscap ?
        capped(input.masses[j],input.filters[j].cap^input.nmax,k/input.nmax) :
        input.masses[j]^(k/input.nmax)) for j in 2:length(mult))
    m=input.masses[1]
    # If the whole first sector lies below lo, the trial moment is below one.
    pred(lo)=capped(m,lo^input.nmax,pwr)+tail(N)<1
    pred(.5) || return nothing
    lo=.5;hi=1.0
    for _ in 1:45
        mid=(lo+hi)/2
        pred(mid) ? (lo=mid) : (hi=mid)
    end
    m>lo^input.nmax || return nothing
    R=tail(N/2)+(m-lo^input.nmax)^(N/(2input.nmax))
    p=sqdef(R/lo^(N/2));Q=first(bounds(input,N));q=1-Q^-2
    (p,q,gap=gap_estimate(p,q,input.beta),trial=N,n0=N÷2,beta=input.beta)
end

function main_diagnostics()
    s8=load_spectra(8);s10=load_spectra(10)
    summaries=Any[]
    for (nmax,mode,eps) in ((8,"ION",0.0),(8,"ION",3e-5),(10,"JN",0.0),
        (10,"JN",3e-5),(10,"ION",1e-3))
        spectra=nmax==8 ? s8 : s10
        best=nothing;minimax=nothing;direct=nothing;counts=0;bestcap=Inf
        for beta in 6.0:.25:16.0
            tab,tails=thermal_table(spectra,beta)
            # Zero-error rows are optimistic mathematical screening, not certified tables.
            maximum(values(tails))<1e-6 || continue
            input=prepare(tab,nmax,eps,mode;beta,allowodd=true)
            bestcap=min(bestcap,input.filters[1].cap)
            for N in 40:4:192
                x=direct_seed(input,N)
                if !isnothing(x) && (isnothing(direct) || x.gap>direct.gap)
                    direct=merge(x,(base_beta=beta,))
                end
                first(bounds(input,N))<2 || continue
                for k in nmax:2:N÷2,r in 1.0:.25:12.0
                    counts+=1
                    row=seed(input,N,k,r)
                    isnothing(row) && continue
                    row=merge(row,(base_beta=beta,))
                    if isnothing(best) || row.gap>best.gap
                        best=row
                    end
                    if isnothing(minimax) || max(row.p,row.q)<max(minimax.p,minimax.q)
                        minimax=row
                    end
                end
            end
        end
        rec=Dict("nmax"=>nmax,"mode"=>mode,"eps"=>eps,"evaluated_lift_candidates"=>counts,
            "minimum_first_cap"=>bestcap,"best_gap"=>best,
            "best_max_defect"=>minimax,"direct"=>direct)
        push!(summaries,rec);println(rec);flush(stdout)
    end
    write_result(joinpath(@__DIR__,"diagnostics.json"),Dict(
        "metadata"=>metadata("adversarial negative-case scan"),"records"=>summaries))
end
main_diagnostics()
