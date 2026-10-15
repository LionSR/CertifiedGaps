include("model.jl")

"Count L=14 symmetry blocks without constructing a Hamiltonian or dense matrix."
function feasibility()
    started = time()
    o = orbit_data(14, 1)
    dims = [length(active_orbits(o, k, 1)) for k in 0:7]
    bytes = 16 * big(maximum(dims))^2
    limit = parse(BigInt, strip(read("/sys/fs/cgroup/memory.max", String)))
    data = Dict("L" => "14", "charge" => "1", "momentum_dimensions" => string.(dims),
        "one_complex_dense_matrix_bytes" => string(bytes),
        "host_memory_limit_bytes" => string(limit), "single_matrix_fits" => bytes <= limit,
        "note" => "Dense L=14 is infeasible here; n0=7 is also odd and ineligible.",
        "wall_seconds" => time()-started,
        "source_hashes" => source_hashes(@__DIR__; dirs = (),
            files = ("model.jl", "feasibility.jl")))
    write_result(joinpath(@__DIR__, "results", "feasibility.json"), data)
end

if abspath(PROGRAM_FILE) == @__FILE__
    feasibility()
end
