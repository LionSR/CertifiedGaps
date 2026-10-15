using Spin2Gap, LinearAlgebra, JSON3, Dates, SHA
const ST = Spin2Gap.Trials
BLAS.set_num_threads(1)
function main()
    source = joinpath(@__DIR__, "../../results/trial-refinement.json")
    input = JSON3.read(read(source, String), Dict{String,Any})
    rows = Any[]
    for row in input["runs"][1:2]
        t = ST.templates(Int.(row["virtual_spins"]))
        x = Float64.(row["x"])
        for batch in 1:5
            x = ST.canonicalize(t, x)
            best = ST.refine(t, x; iterations = 400)
            x = best.x
            println("D=$(t.D) batch=$batch e=$(best.energy) gradient=$(best.gradient_norm)")
            flush(stdout)
        end
        push!(
            rows,
            Dict(
                "virtual_spins"=>t.spins,
                "bond_dimension"=>t.D,
                "parameters"=>length(t.pairs),
                "x"=>x,
                "infinite_mps_energy_per_bond"=>ST.energy(t, x),
            ),
        )
        write_result(
            joinpath(@__DIR__, "../../results/trial-polished.json"),
            Dict(
                "status"=>"numerical-variational-search",
                "certified"=>false,
                "source_input_sha256"=>bytes2hex(sha256(read(source))),
                "runs"=>rows,
                "completed_utc"=>string(now(UTC)),
            ),
        )
    end
end
main()
