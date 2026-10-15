include(joinpath(@__DIR__, "..", "..", "src", "LeanSeeds.jl"))
using .LeanSeeds.Spectrum
using LinearAlgebra, SparseArrays, Random, Printf, SHA, Dates, JSON3

"Count displacement classes for direct masked-shift sparse multiplication."
function packing_counts(n,m)
    states=basis(n,m); lookup=Dict(s=>i for (i,s) in enumerate(states))
    offsets=[Set{Int}() for _ in 1:2n]
    for (i,word) in enumerate(states),j in 1:n
        k=mod1(j+1,n); pj=3^(j-1); pk=3^(k-1)
        sj=mod(div(word,pj),3)-1; sk=mod(div(word,pk),3)-1
        for (t,delta) in enumerate((-1,1))
            abs(sj-delta)<=1 && abs(sk+delta)<=1 || continue
            push!(offsets[2j-2+t],lookup[word+delta*(pk-pj)]-i)
        end
    end
    (;n,m,dimension=length(states),full_dimension=3^n,
      compressed_shift_masks=sum(length,offsets),
      compressed_unique_shifts=length(union(offsets...)),
      full_shift_masks=2n,ratio=3^n/length(states),
      compressed_mask_bytes_64=8length(states)*sum(length,offsets),
      full_mask_bytes_64=8*3^n*2n)
end

rows=[]
for n in 4:12,m in 0:3
    push!(rows,packing_counts(n,m))
end
writecsv(joinpath(DAT,"packing.csv"),rows)
writejson(joinpath(DAT,"packing.json"),rows)
