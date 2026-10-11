using TOML, LinearAlgebra
using Spin2Gap

function main(args)
    length(args) == 2 ||
        error("usage: julia --project=models/spin2 models/spin2/scripts/run.jl CONFIG.toml OUTPUT.json")
    config = TOML.parsefile(args[1])
    BLAS.set_num_threads(get(config, "blas_threads", 1))
    result = Dict{String,Any}("provenance"=>provenance(config), "runs"=>Any[], "complete"=>false)
    write_result(args[2], result)
    mode = config["mode"]
    if mode == "dmrg"
        checkpoint_dir =
            haskey(config, "checkpoint_dir") ?
            normpath(joinpath(dirname(abspath(args[1])), config["checkpoint_dir"])) : nothing
        for n in config["lengths"], chi in config["maxdims"]
            row = dmrg_levels(
                n;
                maxdim = chi,
                nsweeps = config["nsweeps"],
                cutoff = config["cutoff"],
                seed = config["seed"],
                excited = get(config, "excited_sz0", true),
                weight = get(config, "penalty_weight", 10.0),
                outputlevel = get(config, "outputlevel", 1),
                checkpoint_dir,
            )
            if n <= 6
                spectrum = exact_spectrum(n)
                row["full_diagonalization_gap"] = spectrum[2]-spectrum[1]
                row["full_diagonalization_E0"] = spectrum[1]
            end
            push!(result["runs"], row)
            write_result(args[2], result)
        end
    elseif mode == "thermal"
        betas = Float64.(config["betas"])
        all(>(0), betas) || error("seed beta must be positive")
        targets = sort(unique(vcat(betas, 2 .* betas)))
        for n in config["base_lengths"], chi in config["maxdims"], dt in config["dbetas"]
            curves = Dict{Int,Any}()
            for len in (n, 2n)
                curves[len] =
                    thermal_curve(len, targets; dbeta = dt, maxdim = chi, cutoff = config["cutoff"])
            end
            rows = [seed_estimates(curves[n], curves[2n], n, beta) for beta in betas]
            row = Dict(
                "n"=>n,
                "maxdim"=>chi,
                "dbeta"=>dt,
                "seeds"=>rows,
                "curves"=>[
                    Dict("n"=>len, "values"=>[curves[len][b] for b in targets]) for len in (n, 2n)
                ],
            )
            if n <= 6
                es = exact_spectrum(n)
                row["small_chain_max_abs_logZ_error"] =
                    maximum(abs(curves[n][b]["logZ"]-log_partition(es, b)) for b in targets)
            end
            push!(result["runs"], row)
            write_result(args[2], result)
        end
    else
        error("unknown mode $mode")
    end
    result["complete"] = true
    write_result(args[2], result)
    println("Saved numerical results to $(args[2]); certified=false")
end

main(ARGS)
