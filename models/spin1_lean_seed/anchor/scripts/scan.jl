include("../src/common.jl")
using .AnchorStudy, LinearAlgebra, Optim
using CertifiedGaps.Results
const A = AnchorStudy

function run_scan()
    rows = Any[]
    for l in 4:6
        objective(x) = -(2minimum(first(es) for es in values(A.block_spectra(l,x...))) +
            2A.local_min(x...))
        best = nothing
        for x0 in ([0.,0.],[-0.5,0.],[-1.,0.2],[0.,0.5])
            opt = optimize(objective,x0,NelderMead(),Optim.Options(iterations=400,
                x_abstol=1e-7,f_abstol=1e-9))
            if best === nothing || Optim.minimum(opt) < Optim.minimum(best)
                best=opt
            end
        end
        u,v = Optim.minimizer(best)
        # Round the parameters; c is rounded down with a separate safety margin.
        u,v = round(u;digits=4),round(v;digits=4)
        c = floor((A.local_min(u,v)-1e-8)*10^7)/10^7
        E = -objective([u,v])
        println("n=$(2l) optimum u=$u v=$v c=$c H0_ground=$E")
        push!(rows,Dict("n"=>2l,"u"=>u,"v"=>v,"c"=>c,"ground"=>E,
            "bare_ground"=>-objective([0.,0.])))
    end
    write_result(joinpath(A.RESULTS,"optimized.json"),Dict(
        "metadata"=>A.metadata(Dict("optimization"=>"four-start Nelder Mead ground energy")),
        "rows"=>rows))
end
run_scan()
