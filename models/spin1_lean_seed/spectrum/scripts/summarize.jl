include("../src/common.jl")

"""
    filter_runs(records)

One row per filter run and stopping rule, plus per-cycle curves, from the raw filter JSONs; adds
the adaptive schedule (degree 0) that picks the cheapest successful measured run per sector.
"""
function filter_runs(records)
    runs=[]; curves=[]
    for f in records
        for r in f.rows
            push!(curves,(;source=Int(f.source),n=Int(f.n),twist=String(f.twist),
                m=Int(f.m),degree=Int(f.degree),bits=Int(f.bits),k=Int(f.k),p=Int(f.p),
                cycle=Int(r.cycle),total_degree=Int(r.total_degree),hv=Int(r.hv),
                residual=Float64(r.residual),eigerror=Float64(r.eigerror),
                linear_pass=Bool(r.linear_pass),quadratic_pass=Bool(r.quadratic_pass)))
        end
        for method in ("linear","quadratic")
            idx=findfirst(r->getproperty(r,Symbol(method*"_pass")),f.rows)
            r=isnothing(idx) ? last(f.rows) : f.rows[idx]
            push!(runs,(;source=Int(f.source),n=Int(f.n),twist=String(f.twist),
                m=Int(f.m),degree=Int(f.degree),bits=Int(f.bits),method,
                selected_degree=Int(f.degree),
                k=Int(f.k),p=Int(f.p),dimension=Int(f.dimension),
                success=!isnothing(idx),cycles=Int(r.cycle),hv=Int(r.hv),
                residual=Float64(r.residual),eigerror=Float64(r.eigerror),
                target=Float64(getproperty(f.target,Symbol("residual_"*method))),
                seconds=Float64(r.seconds)))
        end
    end
    # Adaptive schedule chooses only measured successful runs, independently per
    # sector. Degree 0 labels the schedule, selected_degree retains the recipe.
    adaptive=[]
    for source in (1,2),bits in (64,128),method in ("linear","quadratic")
        candidates=filter(r->r.source==source && r.bits==bits && r.method==method,runs)
        for key in unique((r.n,r.twist,r.m) for r in candidates)
            options=filter(r->(r.n,r.twist,r.m)==key,candidates)
            successes=filter(r->r.success,options)
            selected=isempty(successes) ? options[argmin([r.hv for r in options])] :
                     successes[argmin([r.hv for r in successes])]
            push!(adaptive,merge(selected,(degree=0,)))
        end
    end
    append!(runs,adaptive)
    sort!(runs;by=r->(r.source,r.bits,r.degree,r.method,r.n,r.twist,r.m))
    runs,curves
end

"""
    summarize()

Rebuild filter_runs.csv, filter_curves.csv, cost_sectors.csv, cost_totals.csv and summary.json.
The raw `filter_s*.json` files are not kept in git; without them the run table is read back from
summary.json, so the cost model (cost_sectors.csv, cost_totals.csv) is still reproducible.
"""
function summarize()
    paths=[p for p in readdir(DAT;join=true)
           if startswith(basename(p),"filter_s") && endswith(p,".json")]
    if isempty(paths)
        runs=[(;(Symbol(k)=>(v isa AbstractString ? String(v) : v) for (k,v) in pairs(r))...)
              for r in JSON3.read(read(joinpath(DAT,"summary.json"),String)).runs]
    else
        runs,curves=filter_runs([JSON3.read(read(p,String)) for p in paths])
        writecsv(joinpath(DAT,"filter_runs.csv"),runs)
        writecsv(joinpath(DAT,"filter_curves.csv"),curves)
    end
    pack=JSON3.read(read(joinpath(DAT,"packing.json"),String))
    costs=[]
    for r in runs
        bits=r.bits==0 ? 64 : r.bits
        factor=r.twist=="C" ? 2 : 1
        logd=ceil(Int,log2(r.dimension))+4
        # Signed dot products need product-sized lanes with accumulation guards.
        gramwidth=2bits+logd
        gram_mults=factor==1 ? 3*r.k*(r.k+1)÷2 : 6*r.k*(r.k+1)
        mixing_entries=r.cycles*r.p^2
        for layout in ("full_ternary","full_Hv_compact_Gram","compressed_direct",
                       "ideal_compressed_24masks")
            D=startswith(layout,"full") ? 3^r.n : r.dimension
            ix=findfirst(x->x.n==r.n && x.m==r.m,pack)
            masks=layout=="compressed_direct" ? Int(pack[ix].compressed_shift_masks) : 2r.n
            hvpasses=factor*(3masks+4ceil(Int,log2(r.n+1)))+(factor==2 ? 12 : 0)
            recurrence_passes=11factor*r.hv
            mixing_passes=4factor*mixing_entries
            # Conditional hybrid: route V and HV once using a fixed permutation
            # network. Factor two covers padding to a power of two lane count.
            compaction_passes=layout=="full_Hv_compact_Gram" ?
                40factor*r.k*ceil(Int,log2(3^r.n)) : 0
            passes=hvpasses*r.hv+recurrence_passes+mixing_passes+compaction_passes
            traffic=BigInt(passes)*bits*D
            gram_bits=gramwidth*(layout=="full_Hv_compact_Gram" ? r.dimension : D)
            # Dimensionless equivalent cost of GMP multiplications at 4M bits.
            mulscale=(gram_bits/4e6)*(log2(gram_bits)/log2(4e6))*
                     (log2(log2(gram_bits))/log2(log2(4e6)))
            gram_bitops=gram_mults*gram_bits*log2(gram_bits)*log2(log2(gram_bits))
            push!(costs,(;source=r.source,n=r.n,twist=r.twist,m=r.m,
                degree=r.degree,bits=r.bits,method=r.method,success=r.success,layout,
                hv=r.hv,vector_bits=bits*D,primitive_passes=passes,
                bit_pass_volume=string(traffic),gram_mults,gram_operand_bits=gram_bits,
                gram_4M_equivalent=gram_mults*mulscale,
                gram_model_bit_operations=gram_bitops,
                total_model_bit_operations=Float64(traffic)+gram_bitops,
                mixing_entries,mixing_recipe_bytes=16factor*mixing_entries,
                small_loop_steps=r.hv*(hvpasses+11factor)+mixing_passes+compaction_passes))
        end
    end
    writecsv(joinpath(DAT,"cost_sectors.csv"),costs)
    totals=[]
    for source in (1,2),bits in (0,64,128),degree in (0,16,32,64),
            method in ("linear","quadratic"),
            layout in ("full_ternary","full_Hv_compact_Gram","compressed_direct",
                       "ideal_compressed_24masks")
        rr=filter(r->r.source==source && r.bits==bits && r.degree==degree &&
                       r.method==method && r.layout==layout,costs)
        isempty(rr) && continue
        push!(totals,(;source,bits,degree,method,layout,sectors=length(rr),
            successes=count(r->r.success,rr),hv=sum(r.hv for r in rr),
            primitive_passes=sum(r.primitive_passes for r in rr),
            bit_pass_volume=string(sum(parse(BigInt,r.bit_pass_volume) for r in rr)),
            gram_mults=sum(r.gram_mults for r in rr),
            gram_4M_equivalent=sum(r.gram_4M_equivalent for r in rr),
            gram_model_bit_operations=sum(r.gram_model_bit_operations for r in rr),
            total_model_bit_operations=sum(r.total_model_bit_operations for r in rr),
            mixing_entries=sum(r.mixing_entries for r in rr),
            mixing_recipe_bytes=sum(r.mixing_recipe_bytes for r in rr),
            small_loop_seconds=sum(r.small_loop_steps for r in rr)/3e4))
    end
    writecsv(joinpath(DAT,"cost_totals.csv"),totals)
    writejson(joinpath(DAT,"summary.json"),merge(metadata(["src/common.jl",
          "src/filtering.jl","scripts/summarize.jl"]),Dict("runs"=>runs,"totals"=>totals)))
end
if abspath(PROGRAM_FILE)==@__FILE__
    summarize()
end
