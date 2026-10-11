using Spin2Gap, LinearAlgebra, JSON3, Dates, SHA
const ST = Spin2Gap.Trials
BLAS.set_num_threads(1)

function main(output = joinpath(@__DIR__, "../results/trial-search.json"))
    rows = Any[]
    for spins in ([1], [1, 2], [1, 2, 3], [0, 1, 2, 3], [0, 1, 2, 3, 4])
        t = ST.templates(spins)
        best = nothing
        started = time()
        for seed in 1:4
            found = ST.search(t; seed, iterations = 3000)
            if isnothing(best) || found.energy < best.energy
                best = found
            end
            println("spins=$spins D=$(t.D) seed=$seed e=$(found.energy) best=$(best.energy)")
            flush(stdout)
        end
        push!(
            rows,
            Dict(
                "virtual_spins"=>spins,
                "bond_dimension"=>t.D,
                "parameters"=>length(t.pairs),
                "pairs"=>[collect(p) for p in t.pairs],
                "x"=>best.x,
                "infinite_mps_energy_per_bond"=>best.energy,
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
main(isempty(ARGS) ? joinpath(@__DIR__, "../results/trial-search.json") : only(ARGS))
