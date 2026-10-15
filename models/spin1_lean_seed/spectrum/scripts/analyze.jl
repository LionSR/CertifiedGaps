include("../src/common.jl")
count_below(v, threshold)=count(x->x<threshold-1e-8,v)

"Infimum cutoff for the dimension-times-exponential tail bound; strict > needed."
function tail_cutoff(values, E0, d, beta, q, relative)
    lo,hi=0.0,4.0
    for _ in 1:70
        c=(lo+hi)/2
        k=count(<(E0+c),values)
        tail=(d-k)*exp(-beta*c)/q
        tail < relative ? (hi=c) : (lo=c)
    end
    hi
end

function analyze()
    sector_rows, eigen_rows, trace_rows, precision_rows = [], [], [], []
    for path in sort(filter(p->startswith(basename(p),"spectrum_"),readdir(DAT;join=true)))
        ring=JSON3.read(read(path,String))
        n,twist,E0=Int(ring.n),String(ring.twist),Float64(ring.E0)
        values=sort(vcat([repeat(Float64.(s.values),s.m==0 ? 1 : 2)
                          for s in ring.sectors]...))
        for s in ring.sectors, m in (s.m==0 ? [0] : [-Int(s.m),Int(s.m)])
            ev=Float64.(s.values)
            push!(sector_rows,(;n,twist,m,dimension=Int(s.dimension),
                lowest=first(ev),computed=length(ev),highest_computed=last(ev),
                k1=count_below(ev,E0+1),k1_5=count_below(ev,E0+1.5),
                k2=count_below(ev,E0+2),k2_5=count_below(ev,E0+2.5),k3=count_below(ev,E0+3),
                residual=Float64(s.residual),method=String(s.method),
                inferred_from_positive_m=m<0))
            for (i,e) in enumerate(ev)
                push!(eigen_rows,(;n,twist,m,index=i,energy=e,gap=e-E0))
            end
        end
        for beta in (12.25,10.5)
            q=sum(exp.(-beta.*(values.-E0)))
            z=q*exp(-beta*(E0+SHIFT*n))
            omitted_z=sum((Int(s.dimension)-length(s.values))*(s.m==0 ? 1 : 2)*
                          exp(-beta*(last(s.values)+SHIFT*n)) for s in ring.sectors)
            cstar=tail_cutoff(values,E0,3^n,beta,q,1e-6)
            push!(trace_rows,(;n,twist,beta,E0,Z_shifted=z,Q_ground_relative=q,
                k1=count_below(values,E0+1),k1_5=count_below(values,E0+1.5),
                k2=count_below(values,E0+2),k2_5=count_below(values,E0+2.5),
                k3=count_below(values,E0+3),cstar,kstar=count(<(E0+cstar),values),
                omitted_z_bound_if_counts_complete=omitted_z,
                tail_c2_relative=(3^n-count(<(E0+2),values))*exp(-2beta)/q,
                tail_c3_relative=(3^n-count(<(E0+3),values))*exp(-3beta)/q))
            for budget in (3e-5,1e-4,1e-3)
                # Width of two-sided enclosure is 2 Z sinh(beta*epsilon).
                eps_all=asinh(budget/(2z))/beta
                eps_half=asinh(budget/(4z))/beta
                cabs=tail_cutoff(values,E0,3^n,beta,q,budget/(4z))
                cutoff=max(cstar,cabs)+1e-8
                for s in ring.sectors
                    ev=Float64.(s.values)
                    k=count(<(E0+cutoff),ev)
                    margin=k==0 ? 0.0 : E0+cutoff-ev[k]
                    push!(precision_rows,(;n,twist,beta,budget,m=Int(s.m),
                        dimension=Int(s.dimension),epsilon_all=eps_all,
                        epsilon_absolute_radius=log1p(budget/z)/beta,
                        epsilon_half=eps_half,cutoff,k,empirical_anchor_margin=margin,
                        residual_linear=eps_half,
                        residual_quadratic=k==0 ? 0.0 : sqrt(eps_half*margin/2)))
                end
            end
        end
    end
    writecsv(joinpath(DAT,"sectors.csv"),sector_rows)
    writecsv(joinpath(DAT,"eigenvalues.csv"),eigen_rows)
    writecsv(joinpath(DAT,"traces.csv"),trace_rows)
    writecsv(joinpath(DAT,"precision.csv"),precision_rows)
    writejson(joinpath(DAT,"analysis.json"),merge(metadata(["src/common.jl","scripts/analyze.jl"]),
        Dict("traces"=>trace_rows,"precision"=>precision_rows)))
    trace_rows,precision_rows
end
if abspath(PROGRAM_FILE)==@__FILE__
    analyze()
end
