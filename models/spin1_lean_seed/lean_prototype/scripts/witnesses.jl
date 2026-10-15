# Write results/witnesses.json: the exact Ritz/residual arithmetic witnesses that the Lean
# kernel checks in src/SeedPrototype/H{6,8,10}.lean, with source hashes and versions.
# This is an inventory, NOT a gap certificate: no physical eigenvalue bound is claimed.
#
#     julia --project=. scripts/witnesses.jl        (after scripts/generate.jl)
using CertifiedGaps.Results: write_result, source_hashes, run_metadata

const MODEL_ROOT = dirname(@__DIR__)

"""Every Lean and Lake source file under src/, relative to the directory root."""
function lean_sources(root=MODEL_ROOT)
    files = String[]
    for (dir, dirs, names) in walkdir(joinpath(root, "src"))
        filter!(!=(".lake"), dirs)
        for name in names
            (endswith(name, ".lean") || name in ("lakefile.toml", "lean-toolchain",
                "lake-manifest.json")) && push!(files, relpath(joinpath(dir, name), root))
        end
    end
    return sort(files)
end

"""Hashes of everything the checker trusts or replays: Lean sources, scripts, tests, check.jl."""
hashed_sources(root=MODEL_ROOT) =
    source_hashes(root; files=vcat(lean_sources(root), ["check.jl"]))

"""Rows of results/exact_witnesses.tsv as exact strings (rationals as `num//den`)."""
function witness_rows(root=MODEL_ROOT)
    lines = readlines(joinpath(root, "results", "exact_witnesses.tsv"))[2:end]
    return map(lines) do line
        f = split(line, '\t')
        Dict("n" => f[1], "degree" => f[2], "lane_bits" => f[3],
            "ritz_upper" => "-$(f[4])//$(f[5])", "residual_upper" => "$(f[6])//$(f[7])",
            "G" => f[8], "A_filter" => f[9], "C_filter" => f[10])
    end
end

function main()
    toolchain = strip(read(joinpath(MODEL_ROOT, "src", "lean-toolchain"), String))
    status = "kernel arithmetic witnesses only; no physical eigenvalue or gap theorem"
    data = Dict(
        "status" => status,
        "witnesses" => witness_rows(),
        "run_metadata" => run_metadata(Dict("lean_toolchain" => toolchain,
            "lean_dependencies" => "Std only (no Mathlib import)");
            source_sha256=hashed_sources(), status, certified=false),
    )
    write_result(joinpath(MODEL_ROOT, "results", "witnesses.json"), data)
end

if abspath(PROGRAM_FILE) == 
    main()
end
