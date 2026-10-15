include("homotopy.jl")
if isempty(ARGS)
    for n in (8,10),g in ("I","P"),sector in 1:3
        continuation(n,g;sector,mode="binned")
    end
else
    for sector in 1:4
        continuation(12,ARGS[1];sector,mode="binned")
    end
end
