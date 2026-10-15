# Numerical artifact audit: fresh exact sector counts, no new spectral certification.
include("src/Suite.jl")
using .LadderSuite
const LS = LadderSuite

function replay_scans()
    start = time()
    dimensions = Dict(L=>LS.sector_dimensions(L) for L in (2,3,4,5,6,8,10))
    records = Any[]
    hashes = Dict{String,String}()
    for tag in ("4p0","3p5","3p0","2p5","2p0","1p5")
        scanfile = joinpath(LS.BASE,"results/scan_J$(tag).json")
        scan = LS.JSON3.read(read(scanfile,String))
        hashes[relpath(scanfile,LS.BASE)] = LS.sha256_file(scanfile)
        logs = Dict{Int,Dict{Rational{Int},Float64}}()
        for L in keys(dimensions)
            path = joinpath(LS.BASE,"results/spectrum_J$(tag)_L$(L).json")
            data = LS.JSON3.read(read(path,String))
            hashes[relpath(path,LS.BASE)] = LS.sha256_file(path)
            actual = Dict((Int(r.m),Int(r.k),Int(r.leg_parity),Int(r.spin_parity))=>
                (Int(r.dimension),Int(r.multiplicity)) for r in data.sectors)
            length(actual) == length(data.sectors) || error("duplicate sector")
            actual == dimensions[L] || error("wrong sector accounting at J=$tag L=$L")
            Int(data.dimension) == 4^L || error("wrong total dimension")
            logs[L] = Dict(parse(Rational{Int},r.beta)=>Float64(r.logZ) for r in data.logZ_grid)
            length(logs[L]) == 128 || error("incomplete trace grid")
        end
        d = (n,beta)->(-expm1(logs[2n][Rational{Int}(beta)]-2logs[n][Rational{Int}(beta)]),
            -expm1(logs[2n][2Rational{Int}(beta)]-2logs[2n][Rational{Int}(beta)]))
        rows = LS.scan_seeds(d,2:5,0.25:0.25:16)
        length(rows) == length(scan.grid) == 256 || error("incomplete scan grid")
        for (fresh,saved) in zip(rows,scan.grid)
            fresh.n₀ == saved.n0 && fresh.β₀ == saved.beta || error("grid mismatch")
            isapprox(fresh.εₛ,saved.eps_s; atol=5e-10,rtol=0) || error("spatial mismatch")
            isapprox(fresh.εₜ,saved.eps_t; atol=5e-10,rtol=0) || error("thermal mismatch")
        end
        best = LS.best_seed(filter(r->iseven(r.n₀),rows))
        @assert (best === nothing) == (scan.best === nothing)
        if best !== nothing
            @assert (best.n₀,best.β₀) == (scan.best.n0,scan.best.beta)
            @assert isapprox(best.bound,scan.best.bound; atol=5e-9,rtol=0)
        end
        push!(records,Dict("J"=>scan.J, "rings"=>sort(collect(keys(dimensions))),
            "largest_dimension"=>"1048576//1", "all_exact_sector_counts_pass"=>true,
            "numerical_grid_reproduced"=>true))
    end
    meta = LS.metadata(Dict("scope"=>"historical numerical scans, fresh exact sector counts"))
    meta["certified"] = false
    meta["status"] = "numerical artifact audit; not a spectral certificate"
    LS.write_result(joinpath(LS.BASE,"results/scan_replay.json"),Dict("metadata"=>meta,
        "artifacts_sha256"=>hashes, "parameters"=>records, "seconds"=>time()-start,
        "limitation"=>"Original dense runs are reused, not rerun; " *
            "trace grids remain floating point."))
    println("PASS: all six historical scans, fresh exact symmetry counts through 20 spins; ",
        time()-start," seconds")
end
replay_scans()
