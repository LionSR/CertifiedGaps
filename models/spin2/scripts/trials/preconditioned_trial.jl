using Spin2Gap, LinearAlgebra, JSON3, SHA, Dates
const SV = Spin2Gap.Variational
BLAS.set_num_threads(1)

function main(
    output = joinpath(@__DIR__, "../../results/trial-preconditioned.json");
    source = joinpath(@__DIR__, "../../results/trial-large.json"),
    multiplicities = nothing,
    batches = 4,
    iterations = 750,
)
    old = last(JSON3.read(read(source, String), Dict{String,Any})["runs"])
    old_spins = Int.(old["virtual_spins"])
    spins =
        isnothing(multiplicities) ? old_spins :
        [j for j in 0:(length(multiplicities)-1) for _ in 1:multiplicities[j+1]]
    model = SV.model(spins)
    x =
        spins == old_spins ? Float64.(old["x"]) :
        SV.embed(
            SV.ST.templates(old_spins),
            Float64.(old["x"]),
            model.t;
            noise = 0.0005,
            seed = 217,
        )
    result = Dict{String,Any}(
        "status"=>"numerical-variational-search",
        "complete"=>false,
        "certified"=>false,
        "uniform_gap_certified"=>false,
        "source_input_sha256"=>bytes2hex(sha256(read(source))),
        "module_sha256"=>bytes2hex(sha256(read(joinpath(@__DIR__, "../../src/Variational.jl")))),
        "runner_sha256"=>bytes2hex(sha256(read(@__FILE__))),
        "runs"=>[
            Dict{String,Any}(
                "virtual_spins"=>model.t.spins,
                "bond_dimension"=>model.t.D,
                "x"=>x,
                "batches"=>Any[],
            ),
        ],
    )
    row = only(result["runs"])
    for batch in 1:batches
        x = SV.ST.canonicalize(model.t, x)
        t0 = time()
        best = SV.optimize_lbfgs(model, x; iterations, tolerance = 1e-8, precondition = true)
        x = best.x
        @assert abs(best.energy-SV.ST.energy(model.t, x)) < 2e-10
        push!(
            row["batches"],
            Dict(
                "batch"=>batch,
                "energy"=>best.energy,
                "gradient_norm"=>best.gradient_norm,
                "reason"=>best.reason,
                "iterations"=>best.iterations,
                "seconds"=>time()-t0,
            ),
        )
        row["x"] = x
        row["energy_per_bond"] = best.energy
        result["complete"] = batch == batches
        result["updated_utc"] = string(now(UTC))
        write_result(output, result)
        println(
            "preconditioned D=$(model.t.D) batch=$batch energy=$(best.energy) gradient=$(best.gradient_norm) $(best.reason) seconds=$(time()-t0)",
        )
        flush(stdout)
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(isempty(ARGS) ? joinpath(@__DIR__, "../../results/trial-preconditioned.json") : only(ARGS))
end
