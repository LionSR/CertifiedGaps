using Spin2Gap, LinearAlgebra, JSON3, SHA, Dates
const SV = Spin2Gap.Variational
BLAS.set_num_threads(1)

function main(output = joinpath(@__DIR__, "../../results/trial-large.json"))
    source = joinpath(@__DIR__, "../../results/trial-analytic.json")
    old = last(JSON3.read(read(source, String), Dict{String,Any})["runs"])
    previous_t = SV.ST.templates(Int.(old["virtual_spins"]))
    previous_x = Float64.(old["x"])
    runs = Any[]
    started = time()
    function save(complete)
        write_result(
            output,
            Dict(
                "status" => "numerical-variational-search",
                "complete" => complete,
                "certified" => false,
                "uniform_gap_certified" => false,
                "runs" => runs,
                "source_input_sha256" => bytes2hex(sha256(read(source))),
                "module_sha256" =>
                    bytes2hex(sha256(read(joinpath(@__DIR__, "../../src/Variational.jl")))),
                "updated_utc" => string(now(UTC)),
                "wall_seconds" => time() - started,
            ),
        )
    end
    for multiplicities in ([4, 6, 7, 6, 3, 1], [5, 8, 9, 8, 5, 2])
        spins = [j for j in 0:(length(multiplicities)-1) for _ in 1:multiplicities[j+1]]
        model = SV.model(spins)
        x = SV.embed(previous_t, previous_x, model.t; noise = 0.0005, seed = 217)
        row = Dict{String,Any}(
            "virtual_spins" => spins,
            "bond_dimension" => model.t.D,
            "parameters" => length(x),
            "multiplicity_dimension" => length(model.t.singlets),
            "equivalent_quartic_terms" => model.equivalent_terms,
            "angular_blocks" => length(model.couplings),
            "batches" => Any[],
            "x" => x,
        )
        push!(runs, row)
        save(false)
        for batch in 1:4
            x = SV.ST.canonicalize(model.t, x)
            t0 = time()
            best = SV.optimize_lbfgs(model, x; iterations = 1500, tolerance = 1e-8)
            x = best.x
            reference = SV.ST.energy(model.t, x)
            @assert abs(reference - best.energy) < 2e-10
            push!(
                row["batches"],
                Dict(
                    "batch" => batch,
                    "energy" => best.energy,
                    "gradient_norm" => best.gradient_norm,
                    "iterations" => best.iterations,
                    "reason" => best.reason,
                    "seconds" => time()-t0,
                    "independent_energy" => reference,
                ),
            )
            row["x"] = x
            row["energy_per_bond"] = best.energy
            save(false)
            println(
                "D=$(model.t.D) batch=$batch energy=$(best.energy) gradient=$(best.gradient_norm) $(best.reason) seconds=$(time()-t0)",
            )
            flush(stdout)
            best.gradient_norm < 1e-8 && break
        end
        previous_t, previous_x = model.t, x
    end
    save(true)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(isempty(ARGS) ? joinpath(@__DIR__, "../../results/trial-large.json") : only(ARGS))
end
