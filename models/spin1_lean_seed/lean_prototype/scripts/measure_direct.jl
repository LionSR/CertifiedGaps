# Sample the actual Lean executable's VmHWM; avoids counting Lake/Julia bootstrap peaks.
function main()
root=dirname(@__DIR__); name=ARGS[1]
cmd=`$(joinpath(root,".tools","lean-4.35.0-rc3-linux","bin","lean")) $(joinpath("SeedPrototype",name*".lean"))`
cmd=addenv(cmd,"LEAN_PATH"=>joinpath(root,".lake","build","lib","lean"))
t0=time_ns(); peak=0; samples=0
p=open(joinpath(root,"benchmarks",name*".direct.txt"),"w") do io
    run(pipeline(ignorestatus(cmd),stdout=io,stderr=io);wait=false)
end
pid=getpid(p)
while process_running(p)
    try
        text=read("/proc/$pid/status",String)
        m=match(r"VmHWM:\s+(\d+) kB",text)
        if m!==nothing
            peak=max(peak,parse(Int,m[1])); samples+=1
        end
    catch e
        isa(e,SystemError) || rethrow()
    end
    sleep(0.02)
end
wait(p)
open(joinpath(root,"benchmarks",name*".direct.time"),"w") do io
    println(io,"wall_seconds=",(time_ns()-t0)/1e9,"\nsampled_lean_peak_kib=",peak,"\nsamples=",samples,"\nexit_code=",p.exitcode)
end
println(name," direct sampled peak KiB=",peak," exit=",p.exitcode)
end
main()
