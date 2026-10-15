# Time the `decide +kernel` benchmark files in lean/SeedKernel/Bench, one fresh Lean process each
# (after `lake build SeedKernel` in lean/).
#
#     julia --startup-file=no scripts/spin1_kernel/benchmark.jl [limit_seconds] [Name ...]
#
# Writes results/spin1_kernel/benchmarks.tsv (all checks) or prints rows (named checks). Raw
# Lean output goes to the git-ignored results/spin1_kernel/bench/. Wall times include Lean
# startup, imports, elaboration and kernel checking; they are not isolated kernel reduction times.
using Printf

const OUT = normpath(joinpath(@__DIR__, "..", "..", "results", "spin1_kernel"))
const LEAN = normpath(joinpath(@__DIR__, "..", "..", "..", "..", "lean"))
const LAKE = get(ENV, "LAKE_BIN", "lake")

"""Benchmark files (in lean/SeedKernel/Bench) with a one-line description, in reporting order."""
const CHECKS = [
    "Baseline" => "imports / trivial check",
    "Bench6" => "n=6 dense matvec + dot",
    "Bench8" => "n=8 dense matvec + dot",
    "Bench10" => "n=10 dense matvec + dot",
    "Bench12" => "n=12 dense matvec + dot",
    "Matvec12" => "n=12 matvec + checksum",
    "Dot12" => "n=12 reversal + dot",
    "Prepared12" => "n=12 prepared masks, matvec + dot",
    "PreparedH8" => "H8 optimized/original moment equality",
    "Gram6" => "n=6 four-vector Gram",
    "PlainDot12" => "n=12 dot, reversed operand supplied",
    "FastDot12" => "n=12 dot with division-free reversal",
    "Fast12" => "n=12 prepared matvec + fast dot",
    "FastReverse8" => "n=8 fast/original reversal equality",
    "PreparedMatvec12" => "n=12 prepared matvec + checksum",
]

"""
    run_one(name, limit) -> (seconds, peak_rss_kib, exit_code)

Run `lake env lean SeedKernel/Bench/<name>.lean` in a fresh child Julia process so that the child
rusage peak belongs to this check alone. Exit code 124 means the time limit was hit.
"""
function run_one(name, limit)
    mkpath(joinpath(OUT, "bench"))
    out = read(`$(Base.julia_cmd()) --startup-file=no $(@__FILE__) --child $name $limit`, String)
    t, rss, code = split(strip(out), '\t')
    return parse(Float64, t), parse(Int, rss), parse(Int, code)
end

function child(name, limit)
    cmd = Cmd(`$LAKE env lean $(joinpath("SeedKernel", "Bench", name * ".lean"))`; dir=LEAN)
    log = joinpath(OUT, "bench", name * ".txt")
    t0 = time_ns()
    p = open(log, "w") do io
        run(pipeline(ignorestatus(cmd); stdout=io, stderr=io); wait=false)
    end
    status = timedwait(() -> !process_running(p), limit)
    status == :timed_out && (kill(p); wait(p))
    elapsed = (time_ns() - t0) / 1e9
    usage = zeros(Clong, 18)
    ccall(:getrusage, Cint, (Cint, Ptr{Clong}), -1, usage) == 0 || error("getrusage failed")
    # ru_maxrss is bytes on macOS and KiB on Linux.
    rss = Sys.isapple() ? usage[5] ÷ 1024 : usage[5]
    code = status == :timed_out ? 124 : p.exitcode
    println(elapsed, '\t', rss, '\t', code)
end

function main(args)
    if !isempty(args) && args[1] == "--child"
        return child(args[2], parse(Float64, args[3]))
    end
    limit = isempty(args) ? 300.0 : parse(Float64, args[1])
    names = length(args) > 1 ? args[2:end] : first.(CHECKS)
    rows = String[]
    for name in names
        t, rss, code = run_one(name, limit)
        @printf("%-17s %8.2f s %6.2f GiB exit %d\n", name, t, rss / 2^20, code)
        push!(rows, join((name, @sprintf("%.2f", t), rss, code), '\t'))
    end
    length(args) > 1 && return
    open(joinpath(OUT, "benchmarks.tsv"), "w") do io
        println(io, "# ", Sys.MACHINE, ", ", Sys.cpu_info()[1].model, ", ",
            strip(read(joinpath(LEAN, "lean-toolchain"), String)), "; exit 124 = over limit")
        println(io, "check\twall_seconds\tpeak_rss_kib\texit_code")
        foreach(r -> println(io, r), rows)
    end
end

main(ARGS)
