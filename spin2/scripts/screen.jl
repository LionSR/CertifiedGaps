using JSON3, SHA
using ThermalMPOGap.Results: write_result

# Necessary beta window if the supplied delta E and gap were exact. For numerical DMRG inputs
# this is only a screening heuristic. Residual norms do not turn these energy differences into
# rigorous interval bounds. Moved verbatim to `ThermalMPOGap.Recurrence.beta_window`.
using ThermalMPOGap.Recurrence: beta_window

function main(args)
    length(args) >= 2 ||
        error("usage: julia --project=spin2 spin2/scripts/screen.jl OUTPUT.json DMRG.json ...")
    rows = Dict{Int,Any}()
    sources = Dict{String,String}()
    for path in args[2:end]
        data = JSON3.read(read(path, String), Dict{String,Any})
        data["complete"] || error("input run is incomplete: $path")
        data["provenance"]["config"]["mode"] == "dmrg" || error("expected DMRG data: $path")
        sources[path] = bytes2hex(sha256(read(path)))
        for row in data["runs"]
            n = row["n"]
            if !haskey(rows, n) || row["maxdim"] > rows[n]["maxdim"]
                rows[n] = row
            end
        end
    end
    screens = Any[]
    for n in sort(collect(keys(rows)))
        haskey(rows, 2n) || continue
        a, b = rows[n], rows[2n]
        delta_energy = b["levels"]["ground_sz0"]["energy"]-2a["levels"]["ground_sz0"]["energy"]
        gaps = filter(!isnothing, [b["sz1_gap_estimate"], b["sz0_gap_estimate"]])
        gap = minimum(gaps)
        screen = beta_window(delta_energy, gap)
        merge!(
            screen,
            Dict(
                "n"=>n,
                "delta_energy_estimate"=>delta_energy,
                "gap_2n_estimate"=>gap,
                "maxdim_n"=>a["maxdim"],
                "maxdim_2n"=>b["maxdim"],
                "ground_residual_n"=>a["levels"]["ground_sz0"]["residual_norm"],
                "ground_residual_2n"=>b["levels"]["ground_sz0"]["residual_norm"],
            ),
        )
        push!(screens, screen)
    end
    isempty(screens) && error("no paired n,2n lengths found")
    output = Dict(
        "status"=>"numerical screening only; no rigorous energy enclosures",
        "certified"=>false,
        "input_sha256"=>sources,
        "screen_script_sha256"=>bytes2hex(sha256(read(@__FILE__))),
        "screens"=>screens,
    )
    write_result(args[1], output)
    for s in screens
        println(
            "n=$(s["n"]): beta >= $(s["beta_min_from_thermal_purity"]), beta <= $(s["beta_max_from_spatial_purity"]); window_nonempty=$(s["window_nonempty"])",
        )
    end
end

abspath(PROGRAM_FILE) == (@__FILE__) && main(ARGS)
