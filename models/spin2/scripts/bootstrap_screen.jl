# The scalar bootstrap is spin-independent. These are targets, not measured spin-2 purities.
using CertifiedGaps.Results: write_result
module BootstrapTools
using CertifiedGaps.ExactArithmetic: frac
using CertifiedGaps.Recurrence: rational_gap_chain
const gap_chain = rational_gap_chain
end
using JSON3, SHA
const BT = BootstrapTools

function main(output = joinpath(@__DIR__, "../results/bootstrap-screen.json"))
    source = joinpath(@__DIR__, "../results/seed-screen.json")
    old = JSON3.read(read(source, String), Dict{String,Any})
    records = Any[]
    for (name, p, q) in (
        ("old_99_percent", BT.frac(1, 100), BT.frac(1, 100)),
        ("balanced_88_percent", BT.frac(3, 25), BT.frac(3, 25)),
        ("asymmetric_95_83_percent", BT.frac(1, 20), BT.frac(17, 100)),
    )
        # n=4 and beta=1 here only fix units for a conditional theorem.
        proof = BT.gap_chain(p, q, 4, BT.frac(1); steps = 6)
        s, t = Float64(1-p), Float64(1-q)
        ceiling = -log(s*t)/log(t/(1-t))
        windows = [
            Dict(
                "n"=>row["n"],
                "beta_min"=>log(t/(1-t))/row["gap_2n_estimate"],
                "beta_max"=>-log(s*t)/row["delta_energy_estimate"],
                "window_nonempty"=>row["delta_energy_over_gap"]<=ceiling,
            ) for row in old["screens"]
        ]
        push!(
            records,
            Dict(
                "target"=>name,
                "spatial_purity_target"=>s,
                "temporal_purity_target"=>t,
                "necessary_deltaE_over_gap_ceiling"=>ceiling,
                "conditional_scalar_proof"=>proof,
                "numerical_windows"=>windows,
            ),
        )
        println(
            "$name: necessary deltaE/gap <= $ceiling; short-chain windows=$(getindex.(windows,"window_nonempty"))",
        )
    end
    write_result(
        output,
        Dict(
            "conditional_scalar_checks_passed"=>true,
            "spin2_purity_seeds_certified"=>false,
            "status"=>"Targets with exact scalar bootstrap checks; DMRG-based windows are numerical screens only",
            "input_sha256"=>bytes2hex(sha256(read(source))),
            "targets"=>records,
        ),
    )
end
abspath(PROGRAM_FILE)==(@__FILE__) &&
    main(isempty(ARGS) ? joinpath(@__DIR__, "../results/bootstrap-screen.json") : only(ARGS))
