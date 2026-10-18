# Consistency replay of the stored floating-point scans in results/scans (NOT a certificate;
# about 1 s). Checks each stored spectrum file's dimension 4^L and complete log Z grid, then
# recomputes all 256 seed rows per coupling from the stored log Z grids. The dense eigensolves
# themselves are not rerun. Writes results/scan_replay.json with the SHA-256 of every replayed
# file.
#
#     julia --project=julia models/heisenberg_ladder/scripts/replay_scans.jl
include(joinpath(@__DIR__, "..", "src", "HeisenbergLadder.jl"))
using .HeisenbergLadder: BASE
using CertifiedGaps.Seeds: scan_seeds, best_seed
using CertifiedGaps.Results: write_result, sha256_file, model_metadata
using JSON3

const TAGS = ("4p0", "3p5", "3p0", "2p5", "2p0", "1p5")
const RINGS = (2, 3, 4, 5, 6, 8, 10)

"Replay all six scans; throws on the first inconsistency."
function replay_scans()
    start = time()
    records = Any[]
    hashes = Dict{String,String}()
    scans = joinpath(BASE, "results", "scans")
    for tag in TAGS
        scanfile = joinpath(scans, "scan_J$(tag).json")
        scan = JSON3.read(read(scanfile, String))
        hashes[relpath(scanfile, BASE)] = sha256_file(scanfile)
        logs = Dict{Int,Dict{Rational{Int},Float64}}()
        for L in RINGS
            path = joinpath(scans, "spectrum_J$(tag)_L$(L).json")
            data = JSON3.read(read(path, String))
            hashes[relpath(path, BASE)] = sha256_file(path)
            Int(data.dimension) == 4^L || error("wrong total dimension at J=$tag L=$L")
            logs[L] = Dict(parse(Rational{Int}, r.beta) => Float64(r.logZ) for r in data.logZ_grid)
            length(logs[L]) == 128 || error("incomplete trace grid")
        end
        d = (n, beta) -> begin
            b = Rational{Int}(beta)
            (-expm1(logs[2n][b] - 2logs[n][b]), -expm1(logs[2n][2b] - 2logs[2n][b]))
        end
        rows = scan_seeds(d, 2:5, 0.25:0.25:16)
        length(rows) == length(scan.grid) == 256 || error("incomplete scan grid")
        for (fresh, saved) in zip(rows, scan.grid)
            fresh.n₀ == saved.n0 && fresh.β₀ == saved.beta0 || error("grid mismatch")
            isapprox(fresh.εₛ, saved.eps_s; atol = 5e-10, rtol = 0) || error("spatial mismatch")
            isapprox(fresh.εₜ, saved.eps_t; atol = 5e-10, rtol = 0) || error("thermal mismatch")
        end
        best = best_seed(filter(r -> iseven(r.n₀), rows))
        (best === nothing) == (scan.best === nothing) || error("best seed presence mismatch")
        if best !== nothing
            (best.n₀, best.β₀) == (scan.best.n0, scan.best.beta0) || error("best seed mismatch")
            isapprox(best.bound, scan.best.bound; atol = 5e-9, rtol = 0) ||
                error("best bound mismatch")
        end
        push!(records, Dict("J" => scan.J, "rings" => collect(RINGS),
            "largest_dimension" => "1048576//1", "numerical_grid_reproduced" => true))
    end
    meta = model_metadata(BASE, Dict("scope" => "numerical scans and their log Z grids"))
    meta["status"] = "numerical artifact audit; not a spectral certificate"
    write_result(joinpath(BASE, "results", "scan_replay.json"), Dict("metadata" => meta,
        "artifacts_sha256" => hashes, "parameters" => records, "seconds" => time() - start,
        "limitation" => "Dense runs are reused, not rerun; trace grids remain floating point."))
    println("PASS: all six stored scans through 20 spins; ", round(time() - start; digits = 1),
        " s")
end

replay_scans()
