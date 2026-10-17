# Why dense L = 14 is out of reach: count the L = 14, Sz = 1 momentum blocks (translation orbits
# whose period is compatible with the momentum) without building any matrix, and compare one
# dense complex block, and the working memory of its verified enclosure, with the host memory
# limit. Floating-point free, but not a certificate.
#
#   julia --project=julia models/spin1_largeD/scripts/feasibility.jl
using CertifiedGaps.Results
include(joinpath(dirname(@__DIR__), "src", "model.jl"))
using .LargeDModel

"Host memory limit: the cgroup limit when one is set (Linux containers), else physical RAM."
function memory_limit()
    path = "/sys/fs/cgroup/memory.max"
    text = isfile(path) ? strip(read(path, String)) : "max"
    return text == "max" ? BigInt(Sys.total_memory()) : parse(BigInt, text)
end

"Momentum block dimensions m = 0..L/2 of the charge-q sector of the ring of L spin-1 sites."
function momentum_dimensions(L, q)
    top = 3^(L - 1)
    shift(c) = (c % top) * 3 + c ÷ top                 # translation by one site
    charge(c) = sum(1 - d for d in digits(c; base = 3, pad = L))   # digit 0, 1, 2 ↦ m = +1, 0, −1
    periods, seen = Int[], falses(3^L)
    for c in 0:(3^L - 1)
        (seen[c+1] || charge(c) != q) && continue
        p, t = 0, c
        while true
            seen[t+1] = true; p += 1; t = shift(t)
            t == c && break
        end
        push!(periods, p)
    end
    return [count(p -> (m * p) % L == 0, periods) for m in 0:(L ÷ 2)]
end

"Write results/feasibility.json with the L = 14 block dimensions and memory comparison."
function feasibility()
    started = time()
    root = dirname(@__DIR__)
    dims = momentum_dimensions(14, 1)
    bytes = 16 * big(maximum(dims))^2
    # verified_eigen on a complex block holds about eight dense N×N complex arrays at its peak
    # (center, eigenvectors, LAPACK workspace, real/imaginary parts, residual and Gram products);
    # the L = 12 runs peaked near that (see README).
    working = 8bytes
    limit = memory_limit()
    data = Dict("L" => "14", "charge" => "1", "momentum_dimensions" => string.(dims),
        "one_complex_dense_matrix_bytes" => string(bytes),
        "host_memory_limit_bytes" => string(limit), "single_matrix_fits" => bytes <= limit,
        "verified_enclosure_working_bytes_estimate" => string(working),
        "verified_enclosure_fits" => working <= limit,
        "note" => "A verified dense L=14 enclosure is infeasible here; n0=7 is also odd and " *
                  "needs a separate transfer-positivity argument.",
        "wall_seconds" => time() - started,
        "metadata" => model_metadata(root, Dict("L" => 14, "charge" => 1); dirs = (),
            files = ("src/model.jl", "scripts/feasibility.jl")))
    write_result(joinpath(root, "results", "feasibility.json"), data)
end

if abspath(PROGRAM_FILE) == @__FILE__
    feasibility()
end
