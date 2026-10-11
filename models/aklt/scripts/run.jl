using TOML, LinearAlgebra
using AKLTPilot

function main(args)
    length(args)==2 || error("usage: julia --project=models/aklt models/aklt/scripts/run.jl CONFIG OUTPUT")
    config = TOML.parsefile(args[1])
    BLAS.set_num_threads(get(config, "blas_threads", 1))
    result = Dict{String,Any}("provenance"=>provenance(config), "runs"=>Any[], "complete"=>false)
    write_result(args[2], result)
    if config["mode"]=="exact"
        graph = honeycomb(config["nx"], config["ny"])
        result["geometry"] = graph
        data = full_spectrum(
            graph.nsites,
            graph.bonds;
            sector_callback = (energies, sectors)->begin
                result["sectors"] = sectors
                result["eigenvalues_counted"] = length(energies)
                write_result(args[2], result)
            end,
        )
        result["spectrum"] = data.energies
        result["ground_energy"] = data.energies[1]
        result["gap"] = data.energies[2]-data.energies[1]
        result["near_zero_count_at_1e-9"] = count(e->abs(e)<1e-9, data.energies)
        result["thermal"] =
            [Dict("beta"=>b, "logZ"=>log_partition(data.energies, b)) for b in config["betas"]]
    elseif config["mode"]=="thermal"
        nx, ny = config["nx"], config["ny"]
        betas = Float64.(config["betas"])
        all(>(0), betas) || error("diagnostic betas must be positive")
        targets = sort(unique(vcat([0.0], betas, 2betas)))
        shapes = get(config, "shapes", [[nx, ny], [2nx, ny], [2nx, 2ny]])
        for chi in config["maxdims"], dt in config["dbetas"]
            row = Dict{String,Any}(
                "maxdim"=>chi,
                "dbeta"=>dt,
                "curves"=>Any[],
                "diagnostics"=>Any[],
                "complete"=>false,
            )
            push!(result["runs"], row)
            curves = Dict{Tuple{Int,Int},Any}()
            for shape in shapes
                x, y = shape
                graph = honeycomb(x, y)
                cr = Dict{String,Any}("geometry"=>graph, "values"=>Any[], "complete"=>false)
                push!(row["curves"], cr)
                curves[(x, y)] = thermal_curve(
                    graph.nsites,
                    graph.bonds,
                    targets;
                    dbeta = dt,
                    maxdim = chi,
                    cutoff = config["cutoff"],
                    move_sites_back_between_gates = get(
                        config,
                        "move_sites_back_between_gates",
                        true,
                    ),
                    callback = values->begin
                        cr["values"] = [values[b] for b in sort(collect(keys(values)))]
                        write_result(args[2], result)
                    end,
                )
                cr["complete"] = true
                write_result(args[2], result)
            end
            if all(haskey(curves, k) for k in ((nx, ny), (2nx, ny), (2nx, 2ny)))
                row["diagnostics"] = [
                    diagnostics(curves[(nx, ny)], curves[(2nx, ny)], curves[(2nx, 2ny)], b) for
                    b in betas
                ]
            end
            row["complete"] = true
            write_result(args[2], result)
        end
    else
        error("unknown mode")
    end
    result["complete"] = true
    write_result(args[2], result)
    println("Saved $(args[2]); certified=false")
end
main(ARGS)
