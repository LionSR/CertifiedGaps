# Timing includes Lean process startup, imports, elaboration and kernel checking.
using Printf
root=dirname(@__DIR__)
name=ARGS[1]; limit=length(ARGS)>1 ? ARGS[2] : "180"
cmd=`timeout $limit $(joinpath(root,".tools","elan","bin","lake")) env lean $(joinpath("SeedPrototype",name*".lean"))`
cmd=addenv(cmd,"ELAN_HOME"=>joinpath(root,".tools","elan"),"PATH"=>joinpath(root,".tools","elan","bin")*":"*ENV["PATH"])
before=zeros(Clong,18)
ccall(:getrusage,Cint,(Cint,Ptr{Clong}),-1,before)
t0=time_ns()
p=open(joinpath(root,"benchmarks",name*".txt"),"w") do io
    run(pipeline(ignorestatus(cmd),stdout=io,stderr=io))
end
elapsed=(time_ns()-t0)/1e9
usage=zeros(Clong,18)
ccall(:getrusage,Cint,(Cint,Ptr{Clong}),-1,usage)==0 || error("getrusage failed")
open(joinpath(root,"benchmarks",name*".time"),"w") do io
    println(io,"wall_seconds=",elapsed,"\nmax_rss_kib=",usage[5],"\npreexisting_child_peak_kib=",before[5],"\nexit_code=",p.exitcode)
end
@printf("%s: %.3fs, %d KiB, exit %d\n",name,elapsed,usage[5],p.exitcode)
exit(p.exitcode)
