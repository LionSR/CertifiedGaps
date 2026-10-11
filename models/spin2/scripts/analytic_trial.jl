using Spin2Gap, LinearAlgebra, JSON3, SHA, Dates
const SV = Spin2Gap.Variational
BLAS.set_num_threads(1)

function main(output = joinpath(@__DIR__, "../results/trial-analytic.json"))
    source = joinpath(@__DIR__, "../results/trial-polished.json")
    old = JSON3.read(read(source, String), Dict{String,Any})["runs"][2]
    previous_t = SV.ST.templates(Int.(old["virtual_spins"]))
    previous_x = Float64.(old["x"])
    runs = Any[]
    started = time()
    models = [
        Int.(old["virtual_spins"]),
        [0, 0, 1, 1, 1, 2, 2, 2, 3, 3],
        [0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 4],
        [0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4],
    ]
    for spins in models
        model = SV.model(spins)
        x =
            spins == previous_t.spins ? previous_x :
            SV.embed(previous_t, previous_x, model.t; noise = 0.005, seed = 207)
        batches = Any[]
        for batch in 1:6
            x = SV.ST.canonicalize(model.t, x)
            t0 = time()
            best = SV.optimize(model, x; iterations = 1500, tolerance = 1e-9)
            x = best.x
            row = Dict(
                "batch"=>batch,
                "energy"=>best.energy,
                "gradient_norm"=>best.gradient_norm,
                "iterations"=>best.iterations,
                "reason"=>best.reason,
                "seconds"=>time()-t0,
            )
            push!(batches, row)
            println(
                "D=$(model.t.D) parameters=$(length(x)) batch=$batch energy=$(best.energy) gradient=$(best.gradient_norm) $(best.reason) seconds=$(row["seconds"])",
            )
            flush(stdout)
            abs(best.energy - SV.ST.energy(model.t, x)) < 2e-10 ||
                error("independent contraction disagreement")
            best.gradient_norm < 1e-9 && break
        end
        previous_t, previous_x = model.t, x
        push!(
            runs,
            Dict(
                "virtual_spins"=>spins,
                "bond_dimension"=>model.t.D,
                "parameters"=>length(x),
                "multiplicity_dimension"=>length(model.t.singlets),
                "equivalent_quartic_terms"=>model.equivalent_terms,
                "x"=>x,
                "batches"=>batches,
                "energy_per_bond"=>first(SV.energy_gradient(model, x)),
            ),
        )
        write_result(
            output,
            Dict(
                "status"=>"numerical-variational-search",
                "certified"=>false,
                "uniform_gap_certified"=>false,
                "runs"=>runs,
                "source_input_sha256"=>bytes2hex(sha256(read(source))),
                "module_sha256"=>bytes2hex(
                    sha256(read(joinpath(@__DIR__, "../src/Variational.jl"))),
                ),
                "completed_utc"=>string(now(UTC)),
                "wall_seconds"=>time()-started,
            ),
        )
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(isempty(ARGS) ? joinpath(@__DIR__, "../results/trial-analytic.json") : only(ARGS))
end
