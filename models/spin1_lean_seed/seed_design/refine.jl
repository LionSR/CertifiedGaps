include("search_new.jl")
function main_refine()
    rows=Any[]
    for eps in (1e-5,3e-5,1e-4,2e-4,3e-4,4e-4,5e-4,7.5e-4,1e-3),
        capN in (false,true), capO in (false,true)
        input=prepare(stored_table(2),10,eps,"ION";allowodd=false,beta=10.5)
        for i in (capN ? () : (3,))
            input.filters[i]=(cap=input.masses[i]^0.1,shift=10,q=[1.0],
                P=vcat(zeros(10),[1.0]),B=input.masses[i])
        end
        if !capO
            input.filters[2]=(cap=input.masses[2]^0.1,shift=10,q=[1.0],
                P=vcat(zeros(10),[1.0]),B=input.masses[2])
        end
        found=scan_seeds(input;trials=(72,96,120,144,168,192),multipliers=1.0:0.125:6.0)
        existing=filter(r->r.trial in (72,120),found)
        rec=Dict("eps"=>eps,"capN"=>capN,"capO"=>capO,
            "filters"=>input.filters,"masses"=>input.masses,
            "best"=>isempty(found) ? nothing : first(found),
            "best_per_trial"=>[first(filter(x->x.trial==n,found))
                for n in unique(x.trial for x in found)],
            "best_existing_trial"=>isempty(existing) ? nothing : first(existing))
        push!(rows,rec)
        println((eps,capN,capO)," best=",isempty(found) ? nothing :
            (first(found).gap,first(found).trial)," existing=",isempty(existing) ? nothing :
            (first(existing).gap,first(existing).p,first(existing).q,
                first(existing).k,first(existing).r))
    end
    write_result(joinpath(@__DIR__,"refinement.json"),Dict(
        "metadata"=>metadata("source 2 error/mass ablation"),"records"=>rows))
end
main_refine()
