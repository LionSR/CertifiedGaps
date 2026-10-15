# Why dense L = 14 is out of reach: count the L = 14, Sz = 1 momentum blocks without building
# any matrix, and compare one dense complex block with the host memory limit.
#
#   julia --project=models/spin1_largeD models/spin1_largeD/scripts/feasibility.jl
include(joinpath(@__DIR__, "..", "model.jl"))

"Host memory limit: the cgroup limit when one is set (Linux containers), else physical RAM."
function memory_limit()
    path = "/sys/fs/cgroup/memory.max"
    text = isfile(path) ? strip(read(path, String)) : "max"
    return text == "max" ? BigInt(Sys.total_memory()) : parse(BigInt, text)
end

"Write results/feasibility.json with the L = 14 block dimensions and memory comparison."
function feasibility()
    started = time()
    root = dirname(@__DIR__)
    o = orbit_data(14, 1)
    dims = [length(active_orbits(o, k, 1)) for k in 0:7]
    bytes = 16 * big(maximum(dims))^2
    limit = memory_limit()
    data = Dict("L" => "14", "charge" => "1", "momentum_dimensions" => string.(dims),
        "one_complex_dense_matrix_bytes" => string(bytes),
        "host_memory_limit_bytes" => string(limit), "single_matrix_fits" => bytes <= limit,
        "note" => "Dense L=14 is infeasible here; n0=7 is also odd and ineligible.",
        "wall_seconds" => time()-started,
        "source_hashes" => source_hashes(root; dirs = (),
            files = ("model.jl", "scripts/feasibility.jl")))
    write_result(joinpath(root, "results", "feasibility.json"), data)
end

if abspath(PROGRAM_FILE) == @__FILE__
    feasibility()
end
