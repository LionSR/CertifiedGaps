# Assemble observed benchmark tables and an explicitly scoped certificate inventory.
using SHA, Printf
const ROOT=dirname(@__DIR__)
function main()
    names=["Baseline","Bench6","Bench8","Bench10","Bench12","H6","H8","H10",
        "Matvec12","Dot12","Prepared12","PreparedH8","Gram6","Tests","PlainDot12","FastDot12","Fast12","FastReverse8","PreparedMatvec12"]
    labels=["Imports / trivial check","n=6 dense matvec + dot","n=8 dense matvec + dot",
        "n=10 dense matvec + dot","n=12 dense matvec + dot","H6, degree 24, 3 moments",
        "H8, degree 32, 3 moments","H10, degree 16, 3 moments","n=12 matvec + checksum",
        "n=12 reversal + dot","n=12 prepared masks, matvec + dot","H8 optimized/original equality",
        "n=6, four-vector Gram","Signed/twisted/floor/Gram cases","n=12 dot, reversed operand supplied","n=12 dot with division-free reversal","n=12 prepared matvec + fast dot","n=8 fast/original reversal equality","n=12 prepared matvec + checksum"]
    lines=["| Check | Wall seconds | Peak RSS GiB | Exit |",
        "|---|---:|---:|---:|"]
    open(joinpath(ROOT,"benchmarks","summary.tsv"),"w") do io
        println(io,"check\twall_seconds\tpeak_rss_kib\tpreexisting_child_peak_kib\texit_code")
        for (name,label) in zip(names,labels)
            fields=Dict(split(line,'=';limit=2)[1]=>split(line,'=';limit=2)[2]
                for line in readlines(joinpath(ROOT,"benchmarks",name*".time")))
            @assert fields["exit_code"]=="0"
            @assert parse(Int,get(fields,"preexisting_child_peak_kib","0")) < 10000
            println(io,join([name,fields["wall_seconds"],fields["max_rss_kib"],
                get(fields,"preexisting_child_peak_kib","0"),fields["exit_code"]],'\t'))
            push!(lines,@sprintf("| %s | %.3f | %.3f | 0 |",label,
                parse(Float64,fields["wall_seconds"]),parse(Int,fields["max_rss_kib"])/2^20))
        end
    end
    report=joinpath(ROOT,"REPORT.md")
    content=read(report,String)
    if occursin("FINAL_BENCHMARK_TABLE",content)
        write(report,replace(content,"FINAL_BENCHMARK_TABLE"=>join(lines,'\n')))
    end
    sources=String[]
    for (dir,dirs,files) in walkdir(ROOT)
        filter!(d -> !(d in (".tools", ".lake", "references", "benchmarks")), dirs)
        rel=relpath(dir,ROOT)
        (startswith(rel,".tools") || startswith(rel,".lake") || startswith(rel,"benchmarks") ||
            startswith(rel,"references")) && continue
        for file in files
            (endswith(file,".lean") || endswith(file,".jl") || endswith(file,".sh") ||
                endswith(file,".toml") || file in ("lean-toolchain","lake-manifest.json")) || continue
            push!(sources,relpath(joinpath(dir,file),ROOT))
        end
    end
    sort!(sources)
    open(joinpath(ROOT,"certificate.json"),"w") do io
        println(io,"{")
        println(io,"  \"status\": \"kernel arithmetic witnesses only; no physical eigenvalue or gap theorem\",")
        println(io,"  \"run_metadata\": {\"julia_version\": \"",VERSION,"\", \"lean_toolchain\": \"",
            strip(read(joinpath(ROOT,"lean-toolchain"),String)),"\", \"mathlib_rev\": \"c55e6e786f49471c72fbddbec5415808896aec1e\"},")
        println(io,"  \"witnesses\": [")
        rows=readlines(joinpath(ROOT,"benchmarks","exact_witnesses.tsv"))[2:end]
        for (i,row) in enumerate(rows)
            f=split(row,'\t')
            println(io,"    {\"n\": \"",f[1],"\", \"degree\": \"",f[2],"\", \"lane_bits\": \"",f[3],
                "\", \"ritz_upper\": \"-",f[4],"//",f[5],"\", \"residual_upper\": \"",f[6],"//",f[7],"\"}",i==length(rows) ? "" : ",")
        end
        println(io,"  ],\n  \"source_hashes\": {")
        for (i,file) in enumerate(sources)
            println(io,"    \"",file,"\": \"",bytes2hex(sha256(read(joinpath(ROOT,file)))),"\"",i==length(sources) ? "" : ",")
        end
        println(io,"  }\n}")
    end
end
main()
