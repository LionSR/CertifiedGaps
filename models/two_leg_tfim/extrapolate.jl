# Exploratory extrapolation only; neither these spectra nor this ansatz are certified.
using CertifiedGaps, CertifiedGaps.Results, CertifiedGaps.Seeds
root=@__DIR__; out=[]
for tag in ("2_1","5_2")
    c=read_exact_json(joinpath(root,"results","scan_g$tag.json"))
    s=read_exact_json(joinpath(root,"results","spectra_g$tag.json"))
    d3=c["summary"][2]["delta"]; d4=c["summary"][3]["delta"]
    xi=1/log(d3/d4)
    ev=s["records"][end]["energies"]; gaps=ev[2:end].-ev[1]
    # δ(n)=δ(4)exp(-(n-4)/xi); dilute excitation weight scales linearly with length.
    grid=collect(1/8:1/8:128)
    w=Dict(beta=>sum(exp.(-beta.*gaps)) for beta in union(grid,2grid))
    function defects(n,beta)
        delta=d4*exp(-(n-4)/xi)
        qn=1+n/8*w[beta]; q2n=1+2n/8*w[beta]
        ps=exp(-beta*delta)*q2n/qn^2
        pt=(1+2n/8*w[2beta])/q2n^2
        (max(0,1-ps),max(0,1-pt))
    end
    first_seed=nothing
    for n in 6:2:128
        b=best_seed(scan_seeds(defects,[n],grid))
        if b!==nothing
            first_seed=b;break
        end
    end
    push!(out,(;g=c["g"],fit_lengths=[3,4],xi_delta=xi,delta4=d4,
        ed_reference_gap=ev[2]-ev[1],first_predicted_even_seed=first_seed,
        caveat="Uncontrolled extrapolation; finite-size gap and dilute-gas ansatz may be biased near criticality."))
end
write_result(joinpath(root,"results","extrapolation.json"),(;status="numerical, not certified",out,
    source_hashes=source_hashes(root;files=["extrapolate.jl"])))
println(out)
