using Spin2Gap, LinearAlgebra, JSON3, Dates, SHA, Random
const ST = Spin2Gap.Trials
BLAS.set_num_threads(1)

function main(output = joinpath(@__DIR__, "../results/trial-refinement.json"))
    initial = JSON3.read(
        read(joinpath(@__DIR__, "../results/trial-search.json"), String),
        Dict{String,Any},
    )
    old = only(filter(r->r["virtual_spins"]==[0, 1, 2, 3], initial["runs"]))
    previous =
        (spins = Int.(old["virtual_spins"]), x = Float64.(old["x"]), pairs = Tuple.(old["pairs"]))
    rows = Any[]
    for spins in ([0, 1, 1, 2, 2, 3], [0, 1, 1, 1, 2, 2, 3], [0, 0, 1, 1, 1, 2, 2, 2, 3, 3])
        t = ST.templates(spins)
        started = time()
        best = nothing
        for seed in 1:3
            rng = MersenneTwister(100+seed)
            x = 0.08*randn(rng, length(t.pairs))
            # Embed every old multiplet in the corresponding new occurrence.
            map = [
                findall(==(j), spins)[count(==(j), previous.spins[1:a])] for
                (a, j) in enumerate(previous.spins)
            ]
            for (i, (a, b)) in enumerate(previous.pairs)
                target = findfirst(==((map[a], map[b])), t.pairs)
                x[target] = previous.x[i]
            end
            found = ST.refine(t, x; iterations = 350)
            if isnothing(best) || found.energy<best.energy
                best=found
            end
            println(
                "spins=$spins D=$(t.D) seed=$seed e=$(found.energy) gradient=$(found.gradient_norm) iterations=$(found.iterations)",
            )
            flush(stdout)
        end
        previous = (spins = spins, x = best.x, pairs = t.pairs)
        push!(
            rows,
            Dict(
                "virtual_spins"=>spins,
                "bond_dimension"=>t.D,
                "parameters"=>length(t.pairs),
                "pairs"=>[collect(p) for p in t.pairs],
                "x"=>best.x,
                "infinite_mps_energy_per_bond"=>best.energy,
                "gradient_norm"=>best.gradient_norm,
                "seconds"=>time()-started,
                "certified"=>false,
            ),
        )
        write_result(
            output,
            Dict(
                "status"=>"numerical-variational-search",
                "runs"=>rows,
                "gap_certified"=>false,
                "completed_utc"=>string(now(UTC)),
                "source_sha256"=>bytes2hex(sha256(read(joinpath(@__DIR__, "../src/Trials.jl")))),
            ),
        )
    end
end
main(isempty(ARGS) ? joinpath(@__DIR__, "../results/trial-refinement.json") : only(ARGS))
