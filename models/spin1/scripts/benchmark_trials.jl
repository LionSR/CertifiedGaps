import Spin1Certificates as SC
import Nemo

function contract_flint(blocks, n)
    records = Any[]
    for block in blocks
        T, O = Nemo.matrix(Nemo.ZZ, block.T), Nemo.matrix(Nemo.ZZ, block.O)
        X = T^(n - 2)
        push!(
            records,
            Dict(
                "delta" => block.delta,
                "dimension" => size(T, 1),
                "U" => string(Nemo.tr(X * O)),
                "V" => string(Nemo.tr(X * T^2)),
            ),
        )
    end
    return sum(SC.integer(r["U"]) for r in records),
    sum(SC.integer(r["V"]) for r in records),
    records
end

labels, A = SC.trial_matrices()
blocks = SC.trial_blocks(labels, A; native = true)
reduced = SC.su2_trial_blocks(labels, A, blocks)
# Warm all three implementations at small N; construction is excluded equally.
SC.contract_trial(blocks, 4)
contract_flint(blocks, 4)
SC.contract_trial_su2(reduced, 4)
outputs = Dict{String,Any}()
for (name, f) in (
    ("julia_bigint_charge", n -> SC.contract_trial(blocks, n)),
    ("flint_charge", n -> contract_flint(blocks, n)),
    ("flint_su2", n -> SC.contract_trial_su2(reduced, n)),
)
    start = time()
    values = [f(n) for n in (72, 120)]
    outputs[name] = (seconds = time() - start, values = values)
    println("$name both lengths: $(outputs[name].seconds)s")
    flush(stdout)
end
@assert outputs["julia_bigint_charge"].values ==
        outputs["flint_charge"].values ==
        outputs["flint_su2"].values
SC.write_result(
    joinpath(@__DIR__, "../results/trial-benchmark.json"),
    Dict(
        "status" => "passed",
        "scope" => "Same process, warmed contraction only, both N=72 and N=120; identical full integers and all charge blocks",
        "julia_bigint_charge_seconds" => outputs["julia_bigint_charge"].seconds,
        "flint_charge_seconds" => outputs["flint_charge"].seconds,
        "flint_su2_seconds" => outputs["flint_su2"].seconds,
        "charge_dimensions" => [size(b.T, 1) for b in blocks],
        "su2_dimensions" => [b.dimension for b in reduced],
        "julia" => string(VERSION),
        "Nemo" => string(pkgversion(Nemo)),
    ),
)
